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
import os

from dataclasses import dataclass

from cadctap import CadcTapClient
from caom2utils import cadc_client_wrapper
from caom2pipe import client_composable as clc
from caom2pipe import manage_composable as mc
from gem2caom2 import obs_metadata as gom
from gem2caom2.obs_file_relationship import GemObsFileRelationship
from gem2caom2.obs_file_relationship import repair_data_label
from gem2caom2 import gem_name, instruments
from gem2caom2.util import Inst, COLLECTION


__all__ = [
    # 'CachingObsFileRelationship',
    'current_instrument',
    'defining_metadata_finder',
    'get_gofr',
    'get_instrument',
    'get_obs_id_from_cadc',
    'get_obs_id_from_headers',
    'get_obs_metadata',
    'init_global',
    'repair_instrument',
    'set_ofr',
]


GEMINI_METADATA_URL = (
    'https://archive.gemini.edu/jsonsummary/canonical/filepre='
)

# lazy initialization for jsonsummary metadata from Gemini
om = None
# lazy initialization for the Gemini listing of files
gofr = None
# value repair cache
value_repair = mc.ValueRepairCache()
# gofr replacement
defining_metadata_finder = None
#
gemini_session = mc.get_endpoint_session()


# globals are BAD, but if this existed in a class instance, that
# class instance would need to be visible in main_app.py, which is not
# something that fits2caom2 handles right now
def get_gofr(config=None):
    global gofr
    if gofr is None:
        gofr = CachingObsFileRelationship(config)
    global defining_metadata_finder
    if defining_metadata_finder is None:
        defining_metadata_finder = DefiningMetadataFinder(config)
    return gofr


def set_ofr(value):
    # replace the lookup functionality of the file listing from Gemini,
    # with the query results from archive.gemini.edu
    global gofr
    gofr = value


def init_global(config):
    global om
    om = gom.json_lookup
    get_gofr(config)


def get_obs_metadata(file_id):
    """
    Download the Gemini observation metadata for the given obs_id.

    :param file_id: The file ID
    :return: Dictionary of observation metadata.
    """
    logging.debug(f'Begin get_obs_metadata for {file_id}')
    global om
    if om.contains(file_id):
        om.reset_index(file_id)
    else:
        # for TaskType.SCRAPE
        if gofr.query_session is None:
            logging.warning(f'No external access. No observation metadata.')
        else:
            gemini_url = f'{GEMINI_METADATA_URL}{file_id}'

            # Open the URL and fetch the JSON document for the observation
            response = None
            try:
                response = mc.query_endpoint_session(
                    gemini_url, gemini_session
                )
                metadata = response.json()
            finally:
                if response is not None:
                    response.close()
            if len(metadata) == 0:
                raise mc.CadcException(
                    f'Could not find JSON record for {file_id} at '
                    f'archive.gemini.edu.'
                )
            om.add(metadata, file_id)
    logging.debug(f'End get_obs_metadata for {file_id}')


class CachingObsFileRelationship(GemObsFileRelationship):
    """
    The locations in which the relationship between an observationID (data
    label) and a file name can be determined are:
        'connected' order, no files on disk (i.e. not doing TaskType.SCRAPE):
        1 - the specifically constructed file from Paul, which is time-limited
        2 - CAOM entries
        3 - archive.gemini.edu entries

        'connected' order, files on disk (i.e. not doing TaskType.SCRAPE):
        1 - the specifically constructed file from Paul, which is time-limited
            (since it's already in memory, it should be faster to access)
        2 - file on disk
        3 - CAOM entries
        4 - archive.gemini.edu entries

        'unconnected order':
        1 - Paul's file
        2 - file on disk

    This class queries each of these locations in the declared order, and
    then if the answer comes from either 2 or 3, adds to the cache of the
    file.
    """

    def __init__(self, config=None):
        super(CachingObsFileRelationship, self).__init__()
        self._use_local_files = False
        self._is_connected = True
        self._collection = 'GEMINI'
        if config is not None:
            self._use_local_files = config.use_local_files
            self._is_connected = config.is_connected
            self._collection = config.collection
        # use accessor methods for _tap_client, because of how this class
        # will eventually be used - as a global, accessible by all and
        # everywhere, and initialized before there's a config
        self._tap_client = None
        self._session = None
        if self._is_connected:
            self._session = mc.get_endpoint_session()
        self._logger = logging.getLogger(__name__)

    @property
    def query_session(self):
        return self._session

    @property
    def tap_client(self):
        return self._tap_client

    @tap_client.setter
    def tap_client(self, value):
        self._tap_client = value


def get_obs_id_from_cadc(
    file_id, tap_client, collection='GEMINI', update_cache=None
):
    logging.debug(f'Begin get_obs_id_from_cadc for {file_id}')
    file_name = gem_name.GemName.get_file_name_from(file_id)
    query_string = f"""
    SELECT O.observationID, A.lastModified
    FROM caom2.Observation AS O
    JOIN caom2.Plane AS P on P.obsID = O.obsID
    JOIN caom2.Artifact AS A on A.planeID = P.planeID
    WHERE A.uri LIKE '%{file_name}'
    AND O.collection = '{collection}'
    """
    table = clc.query_tap_client(query_string, tap_client)
    result = None
    if len(table) == 1:
        result = table[0]['observationID']
        ut_datetime_str = table[0]['lastModified']
        if update_cache is not None:
            update_cache(file_id, result, ut_datetime_str)
    logging.debug(f'End get_obs_id_from_cadc {result}')
    return result


def get_obs_id_from_headers(file_id, update_cache=None):
    logging.debug(f'Begin get_obs_id_from_headers for {file_id}')
    temp = gem_name.GemName.remove_extensions(file_id)
    try_these = [
        f'{os.getcwd()}/{temp}.fits',
        f'{os.getcwd()}/{temp}.fits.header',
        f'{os.getcwd()}/{temp}.fits.bz2',
        f'{os.getcwd()}/{temp}.fits.gz',
    ]
    result = None
    for f_name in try_these:
        if os.path.exists(f_name):
            headers = cadc_client_wrapper.get_local_file_headers(f_name)
            result = headers[0].get('DATALAB').upper()
            if result is not None:
                if update_cache is not None:
                    update_cache(file_id, result, headers[0].get('DATE'))
                break
    logging.debug(f'End get_obs_id_from_headers {result}')
    return result


@dataclass
class DefiningMetadata:
    """The metadata that is required to know 'what to do next' for
    ingesting a file.

    The instrument changes the rules by which the observation ID is set, and
    thus changes the cardinality handling. Historical Note: It is the need
    for the value of the instrument going forward that makes the large file
    containing file names, data labels, checksums, dates, no longer useful
    for ingestion, and so, the use of this file (/app/data/from_paul.txt)
    was removed from the pipeline.
    """
    instrument: Inst
    data_label: str


class DefiningMetadataFinder:
    """
    A class to handle the hierarchy of querying for metadata, with
    the goal of the least amount of load on archive.gemini.edu.

    The logic implemented in get:

      Looking for provenance | Connected | Use Local
    0          F                   F           F
    1          F                   F           T
    2          F                   T           F
    3          F                   T           T
    4          T                   F           F
    5          T                   F           T
    6          T                   T           F
    7          T                   T           T

    If Connected == F, may only check local (cases 0, 1, 4, 5).
    If Use Local == T, check local (3, 7)
    If Use Local == F, check cadc first, then check gemini (2, 6)

    If Use Local == T, it does no harm to also check cadc, then gemini for
    metadata, because the likelihood is the metadata search will not go
    beyond this point in the execution. This behaviour covers the Looking
    for provenance T/F case.
    """

    def __init__(self, config):
        self._use_local_files = config.use_local_files
        self._connected = mc.TaskType.SCRAPE not in config.task_types
        self._json_lookup = gom.json_lookup
        subject = clc.define_subject(config)
        self._tap_client = CadcTapClient(
            subject=subject, resource_id=config.tap_id
        )
        self._logger = logging.getLogger(self.__class__.__name__)

    def _check_caom2(self, uri):
        self._logger.debug(f'Begin _check_caom2 for {uri}')
        query_string = f"""
        SELECT O.observationID, O.instrument_name
        FROM caom2.Observation AS O
        JOIN caom2.Plane AS P on P.obsID = O.obsID
        JOIN caom2.Artifact AS A on A.planeID = P.planeID
        WHERE A.uri = '{uri}' 
        AND O.collection = '{COLLECTION}'
        """
        table = clc.query_tap_client(query_string, self._tap_client)
        result = None
        if len(table) == 1:
            result = DefiningMetadata(
                Inst(table[0]['instrument_name']), table[0]['observationID'],
            )
        self._logger.debug('End _check_caom2')
        return result

    def _check_local(self, f_name):
        self._logger.debug(f'Begin _check_local for {f_name}')
        file_id = gem_name.GemName.remove_extensions(f_name)
        try_these = [
            f'{os.getcwd()}/{file_id}.fits',
            f'{os.getcwd()}/{file_id}.fits.header',
            f'{os.getcwd()}/{file_id}.fits.bz2',
            f'{os.getcwd()}/{file_id}.fits.gz',
        ]
        result = None
        for f_name in try_these:
            if os.path.exists(f_name):
                headers = cadc_client_wrapper.get_local_file_headers(f_name)
                temp = headers[0].get('DATALAB').upper()
                if temp is not None:
                    result = DefiningMetadata(
                        repair_instrument(headers[0].get('INSTRUME')),
                        headers[0].get('DATALAB')
                    )
                    break
        self._logger.debug('End _check_local')
        return result

    def _check_remote(self, f_name):
        self._logger.debug(f'Begin _check_remote for {f_name}')
        # using the global om structure to look up and store
        # metadata will modify the internal index of the class - maintain
        # that index here with a save/restore, as the lookup can occur for
        # provenance metadata
        looking_for_file_id = gem_name.GemName.remove_extensions(f_name)
        # first check the cache
        current_file_id = self._json_lookup.current
        if self._json_lookup.contains(looking_for_file_id):
            self._json_lookup.reset_index(looking_for_file_id)
        else:
            # this is a remote call to archive.gemini.edu, which updates
            # the json_lookup cache of metadata.
            get_obs_metadata(looking_for_file_id)
        result = None
        if self._json_lookup.contains(looking_for_file_id):
            result = DefiningMetadata(
                Inst(repair_instrument(self._json_lookup.get('instrument'))),
                self._json_lookup.get('data_label'),
            )
        if current_file_id is not None:
            self._json_lookup.reset_index(current_file_id)
        self._logger.debug('End _check_remote')
        return result

    def get(self, uri):
        """

        :param uri: Artifact URI at CADC
        :return: DefiningMetadata instance
        """
        ignore_scheme, ignore_collection, f_name = mc.decompose_uri(uri)
        if self._connected:
            result = None
            if self._use_local_files:
                result = self._check_local(f_name)
            if result is None:
                result = self._check_caom2(uri)
            if result is None:
                result = self._check_remote(f_name)
        else:
            result = self._check_local(f_name)
        if result is not None:
            repaired_data_label = repair_data_label(f_name, result.data_label)
            result.data_label = repaired_data_label
        return result


def repair_instrument(in_name):
    if in_name == 'ALOPEKE':
        # because the value in JSON is a different case than the value in
        # the FITS header
        in_name = 'Alopeke'
    if in_name == 'ZORRO':
        in_name = 'Zorro'
    return Inst(in_name)


# a globally accessible pointer to the latest instance of InstrumentType,
# placed here so that it survives the importlib.import_module call from
# fits2caom2
current_instrument = None


def get_instrument(uri):
    dm = defining_metadata_finder.get(uri)
    global current_instrument
    current_instrument = instruments.instrument_factory(dm.instrument)
    return dm.instrument
