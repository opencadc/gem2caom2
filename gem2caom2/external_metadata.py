# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#

import logging
import re
import requests
from enum import Enum
from requests.adapters import HTTPAdapter
from urllib3 import Retry

from bs4 import BeautifulSoup

from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from gem2caom2.svofps import filter_metadata
from gem2caom2 import gemini_obs_metadata as gom
from gem2caom2.obs_file_relationship import GemObsFileRelationship
from gem2caom2 import gem_name


__all__ = ['get_gofr', 'Inst', 'get_obs_metadata', 'get_pi_metadata',
           'get_filter_metadata', 'set_ofr', 'init_global',
           'CachingObsFileRelationship']


GEMINI_METADATA_URL = \
    'https://archive.gemini.edu/jsonsummary/canonical/filepre='

# lazy initialization for jsonsummary metadata from Gemini
om = None
# lazy initialization for filter metadata from SVO
fm = {}
# lazy initialization for program metadata from Gemini
pm = {}
# lazy initialization for the Gemini listing of files
gofr = None


def get_gofr():
    global gofr
    if gofr is None:
        gofr = GemObsFileRelationship()
    return gofr


def set_ofr(value):
    # replace the lookup functionality of the listing of file from Gemini,
    # with the query results from archive.gemini.edu
    global gofr
    gofr = value


def init_global(incremental):
    global om
    if incremental:
        om = gom.GeminiObsMetadataIncremental()
    else:
        om = gom.GeminiObsMetadata()


class Inst(Enum):

    BHROS = 'bHROS'
    CIRPASS = 'CIRPASS'
    F2 = 'F2'
    FLAMINGOS = 'FLAMINGOS'
    GMOS = 'GMOS'
    GMOSN = 'GMOS-N'
    GMOSS = 'GMOS-S'
    GNIRS = 'GNIRS'
    GPI = 'GPI'
    GRACES = 'GRACES'
    GSAOI = 'GSAOI'
    HOKUPAA = 'Hokupaa+QUIRC'
    HRWFS = 'hrwfs'
    MICHELLE = 'michelle'
    NICI = 'NICI'
    NIFS = 'NIFS'
    NIRI = 'NIRI'
    OSCIR = 'OSCIR'
    PHOENIX = 'PHOENIX'
    TEXES = 'TEXES'
    TRECS = 'TReCS'


def get_obs_metadata(file_id):
    """
    Download the Gemini observation metadata for the given obs_id.

    :param file_id: The file ID
    :return: Dictionary of observation metadata.
    """
    logging.debug('Begin get_obs_metadata for {}'.format(file_id))
    global om
    if om.contains(file_id):
        om.reset_index(file_id)
    else:
        gemini_url = '{}{}'.format(GEMINI_METADATA_URL, file_id)

        # Open the URL and fetch the JSON document for the observation
        session = requests.Session()
        retries = 10
        retry = Retry(total=retries, read=retries, connect=retries,
                      backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        try:
            response = session.get(gemini_url, timeout=20)
            metadata = response.json()
            response.close()
        except Exception as e:
            raise mc.CadcException(
                f'Unable to download Gemini observation metadata from '
                f'{gemini_url} because {str(e)}')
        om.add(metadata, file_id)
    logging.debug('End get_obs_metadata for {}'.format(file_id))


def get_pi_metadata(program_id):
    global pm
    if program_id in pm:
        metadata = pm[program_id]
    else:
        program_url = 'https://archive.gemini.edu/programinfo/' + program_id

        # Open the URL and fetch the JSON document for the observation
        session = requests.Session()
        retries = 10
        retry = Retry(total=retries, read=retries, connect=retries,
                      backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        try:
            response = session.get(program_url, timeout=20)
            xml_metadata = response.text
            response.close()
        except Exception as e:
            raise mc.CadcException(
                'Unable to download Gemini observation metadata from {} '
                'because {}'.format(program_url, str(e)))
        metadata = None
        soup = BeautifulSoup(xml_metadata, 'lxml')
        tds = soup.find_all('td')
        if len(tds) > 0:
            # sometimes the program id points to an html page with an empty
            # table, see e.g. N20200210S0077_bias
            title = None
            if len(tds[1].contents) > 0:
                title = tds[1].contents[0].replace('\n', ' ')
            pi_name = None
            if len(tds[3].contents) > 0:
                pi_name = tds[3].contents[0]
            metadata = {'title': title,
                        'pi_name': pi_name}
            pm[program_id] = metadata
        logging.debug('End get_obs_metadata')
    return metadata


def get_filter_metadata(instrument, filter_name):
    """A way to lazily initialize all the filter metadata reads from SVO."""
    global fm
    repaired_inst = _repair_instrument_name_for_svo(instrument)
    repaired_filters = _repair_filter_name_for_svo(instrument, filter_name)
    if repaired_filters is None:
        # nothing to look up, try something else
        return None
    if repaired_inst in fm and repaired_filters in fm[repaired_inst]:
        result = fm[repaired_inst][repaired_filters]
        if result is not None:
            result.adjust_resolving_power()
    else:
        result = filter_metadata(repaired_inst, repaired_filters)
        if repaired_inst in fm:
            temp = fm[repaired_inst]
            temp[repaired_filters] = result
        else:
            fm[repaired_inst] = {repaired_filters: result}
    return result


def _repair_instrument_name_for_svo(instrument):
    """
    Instrument names from JSON/headers are not necessarily the same
    as the instrument names used by the SVO Filter service. Correlate
    the two here.
    :param instrument the Gemini version
    :return instrument the SVO version
    """
    result = instrument.value
    if instrument is Inst.HRWFS:
        telescope = om.get('telescope')
        if telescope is None:
            obs_id = om.get('data_label')
            raise mc.CadcException(
                '{}: No observatory information for {}'.format(instrument,
                                                               obs_id))
        else:
            if 'Gemini-South' == telescope:
                result = 'AcqCam-S'
            else:
                result = 'AcqCam-N'
    elif instrument is Inst.F2:
        result = 'Flamingos2'
    elif instrument is Inst.FLAMINGOS:
        result = 'Flamingos'
    return result


def _repair_filter_name_for_svo(instrument, filter_names):
    """
    Filter names from JSON/headers are not necessarily the same
    as the filter names used by the SVO Filter service. Correlate
    the two here.

    DB - 02-04-19 - strip the bar code from the filter names

    :param instrument what repairs to apply
    :param filter_names the Gemini version, which may include multiple names
        separated by '+'
    :return filter_name the SVO version
    """
    FILTER_REPAIR_NICI = {'CH4-H4S': 'ED451',
                          'CH4-H4L': 'ED449',
                          'CH4-H1S': 'ED286',
                          'CH4-H1Sp': 'ED379',
                          '': 'ED299',
                          'CH4-H1L': 'ED381',
                          'CH4-H1L_2': 'ED283'}
    # note the lookup repair values are not what comes from the files,
    # they're what's left after the re.sub calls have completed
    # DB 06-05-19
    # The NIRI filter should map to SVO’s NIRI.CO2-0bh-G0225. bh = band-head.
    FILTER_REPAIR_NIRI = {'H2v=2-1s1-G0220': 'H2S1v2-1-G0220',
                          'H2v=2-1S1-G0220w': 'H2S1v2-1-G0220w',
                          'H2v=2-1S1-G0220': 'H2S1v2-1-G0220',
                          'H2v=1-0s1-G0216': 'H2S1v1-0-G0216',
                          'H2v=1-0S1-G0216': 'H2S1v1-0-G0216',
                          'H2Oice_G0230': 'H2Oice-G0230w',
                          'Brgamma-G0218': 'BrG-G0218',
                          'Bra-G0238': 'BrAlpha-G0238',
                          'Bracontt-G0237': 'BrAlphaCont-G0237',
                          'CH4ice227-G0243': 'CH4ice2275-G0243',
                          'COv=2-0bh-G0225': 'CO2-0bh-G0225',
                          'hydrocarb-G0231': 'hydrocarbon-G0231',
                          'H2Oice204-G0242': 'H2Oice2045-G0242',
                          'Jcont121-G0232': 'Jcont1207-G0232',
                          'H2v=2-1s1_G0220': 'H2S1v2-1-G0220'}
    # DB 23-04-19
    # The Qs-18.3um is likely intended to be the same as Qa since 18.3 is the
    # central wavelength of that filter.
    FILTER_REPAIR_TRECS = {'K': 'k',
                           'L': 'l',
                           'M': 'm',
                           'N': 'n',
                           'Nprime': 'nprime',
                           'Qw': 'Qwide',
                           'Qs': 'Qa',
                           'NeII_ref2': 'NeII_ref'}
    FILTER_REPAIR_MICHELLE = {'I79B10': 'Si1',
                              'I88B10': 'Si2',
                              'I97B10': 'Si3',
                              'I103B10': 'Si4',
                              'I105B53': 'N',
                              'I112B21': 'Np',
                              'I116B9': 'Si5',
                              'I125B9': 'Si6',
                              'I185B9': 'Qa',
                              'I209B42': 'Q'}
    # DB 02-04-19
    # The GSAOI filter CO2360 should map to SVO filter GSAOI.CO
    # DB 04-24-19
    # H2(1-0) filter maps to SVO H2_1-0
    # DB 04-30-19
    # Kcntshrt and HeI-2p2s should map to Kshort_cont and HeI2p2s.
    # PaG for GS-CAL20180731-5-017 = SVO’s HIPaGamma
    # DB 02-05-19
    # BrG is Brackett Gamma again, so in SVO it is GSAOI.HIBrGamma.
    # I think H2(2-1) must be GSAOI.H2_2-1_S1.  That observation shows
    # up when you search the Gemini archive for that particular
    # observation and set the filter to “H2 2-1 (S1)“.
    # DB 06-05-19
    # PaB = HIPaBeta.
    FILTER_REPAIR_GSAOI = {'BrG': 'HIBrGamma',
                           'CO2360': 'CO',
                           'HeI1083': 'HeI',
                           'HeI-2p2s': 'HeI2p2s',
                           'H2(1-0)': 'H2_1-0',
                           'H2(2-1)': 'H2_2-1_S1',
                           'Kcntlong': 'Klong_cont',
                           'Kcntshrt': 'Kshort_cont',
                           'PaB': 'HIPaBeta',
                           'PaG': 'HIPaGamma'}

    result = []
    for filter_name in filter_names.split('+'):
        temp = filter_name
        if instrument is Inst.NIRI:
            temp = re.sub(r'con', 'cont', temp)
            temp = re.sub(r'_', '-', temp)
            temp = re.sub('\\(', '', temp)
            temp = re.sub('\\)', '', temp)
            if temp in FILTER_REPAIR_NIRI:
                temp = FILTER_REPAIR_NIRI[temp]
        elif instrument is Inst.NICI:
            if temp in FILTER_REPAIR_NICI:
                temp = FILTER_REPAIR_NICI[temp]
            else:
                logging.info(
                    '{} filter {} not at SVO.'.format(instrument, temp))
                temp = None
        elif instrument is Inst.TRECS:
            temp = filter_name.split('-')
            if len(temp) > 0:
                temp = temp[0]
            if temp in FILTER_REPAIR_TRECS:
                temp = FILTER_REPAIR_TRECS[temp]
        elif instrument is Inst.MICHELLE:
            temp = filter_name.split('-')
            if len(temp) > 0:
                temp = temp[0]
            if temp in FILTER_REPAIR_MICHELLE:
                temp = FILTER_REPAIR_MICHELLE[temp]
        elif instrument is Inst.HRWFS:
            # “ND” in the filter name means ‘neutral density’.  Ignore any
            # of these as they have no impact on the transmitted wavelengths
            # - I think #159 was the only one delivered according to
            # http://www.gemini.edu/sciops/telescope/acqcam/acqFilterList.html.
            # Acqcam/hrwfs was used mainly to look for rapid variability in
            # bright, stellar objects that were really too bright for an 8'
            # telescope and would have saturated the detector without an ND
            # filter.
            if temp.startswith('ND'):
                continue
            temp = temp[0]
        elif instrument is Inst.GSAOI:
            if temp in FILTER_REPAIR_GSAOI:
                temp = FILTER_REPAIR_GSAOI[temp]
        elif instrument is Inst.F2:
            if temp == 'J-lo':
                temp = 'Jlow'
        if temp is not None:
            result.append(temp)
    if len(result) > 0:
        return '+'.join(i for i in result)
    else:
        return None


class CachingObsFileRelationship(GemObsFileRelationship):
    """
    The locations in which the relationship between an observationID (data
    label) and a file name can be determined are:
    1 - the specifically constructed file from Paul, which is time-limited
    2 - CAOM entries
    3 - archive.gemini.edu entries

    This class queries each of these locations in the declared order, and
    then if the answer comes from either 2 or 3, adds to the cache of the
    file.
    """

    def __init__(self):
        super(CachingObsFileRelationship, self).__init__()
        # use accessor methods for _tap_client, because of how this class
        # will eventually be used - as a global, accessible by all and
        # everywhere, and initialized before there's a config
        self._tap_client = None
        self._logger = logging.getLogger(__name__)

    @property
    def tap_client(self):
        return self._tap_client

    @tap_client.setter
    def tap_client(self, value):
        self._tap_client = value

    def get_obs_id(self, file_name):
        file_id = gem_name.GemName.remove_extensions(file_name)
        result = super(CachingObsFileRelationship, self).get_obs_id(file_id)
        if result is None:
            result = self._get_obs_id_from_cadc(file_name, file_id)
            if result is None:
                result = self._get_obs_id_from_gemini(file_id)
        return result

    def _get_obs_id_from_cadc(self, file_name, file_id):
        self._logger.debug(f'Begin _get_obs_id_from_cadc for {file_id}')
        artifact_uri = cc.build_artifact_uri(
            file_name, gem_name.COLLECTION, gem_name.SCHEME)
        query_string = f"""
        SELECT O.observationID 
        FROM caom2.Observation AS O
        JOIN caom2.Plane AS P on P.obsID = O.obsID
        JOIN caom2.Artifact AS A on A.planeID = P.planeID
        WHERE A.artifact_uri = '{artifact_uri}'
        """
        table = mc.query_tap_client(query_string, self._tap_client)
        result = None
        if len(table) == 1:
            result = table[0]['observationID']
            # TODO don't know if the timestamp needs to be useful
            self._update_cache(file_id, result, 1.0)
        self._logger.debug('End _get_obs_id_from_cadc')
        return result

    def _get_obs_id_from_gemini(self, file_id):
        # using the global om structure to look up and store
        # metadata will modify the internal index of the class - maintain
        # that index here with a save/restore
        self._logger.debug(f'Begin _get_obs_id_from_gemini for {file_id}')
        global om
        current_file_id = om.current
        get_obs_metadata(file_id)
        obs_id = om.get('data_label')
        ut_datetime_str = om.get('lastmod')
        ut_datetime_s = mc.make_seconds(ut_datetime_str)
        self._update_cache(file_id, obs_id, ut_datetime_s)
        om.reset_index(current_file_id)
        self._logger.error(f'End _get_obs_id_from_gemini.')
        return obs_id

    def _update_cache(self, file_id, obs_id, dt_s):
        # array contents are:
        # 0 - data label / observation ID
        # 1 - timestamp
        mc.append_as_array(self.name_list, file_id, [obs_id, dt_s])
