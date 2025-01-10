# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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
#  : 4 $
#
# ***********************************************************************
#

"""
archive.gemini.edu provides a JSON record for every archive record. As that
JSON record is the preferred source of metadata (it's usually most correct
than the FITS headers), any ingestion will require retrieving that record.

This module is the classes and methods that do and use the retrieval.
"""


import logging
from astropy import units
from astropy.time import TimeDelta
from collections import deque
from os import path

from cadcutils import exceptions
from cadcdata import FileInfo
from caom2utils import data_util
from caom2pipe.astro_composable import get_datetime_mjd
from caom2pipe import client_composable as clc
from caom2pipe.execute_composable import MetaVisitRunnerMeta, OrganizeExecutesRunnerMeta
from caom2pipe import manage_composable as mc
from gem2caom2.util import Inst
from gem2caom2 import obs_file_relationship


__all__ = [
    'GeminiMetadataLookup',
    'GEMINI_METADATA_URL',
    'HEADER_URL',
    'ProvenanceFinder',
    'repair_instrument',
    'retrieve_json',
]


GEMINI_METADATA_URL = 'https://archive.gemini.edu/jsonsummary/canonical/filepre='
HEADER_URL = 'https://archive.gemini.edu/fullheader/'


class GeminiMetadataLookup:
    def __init__(self, metadata_reader):
        self._reader = metadata_reader
        self._max_exputend = {}

    def camera(self, uri):
        return self._search_json(uri, 'camera')

    def central_wavelength(self, uri):
        return self._search_json(uri, 'central_wavelength')

    def data_label(self, uri):
        instrument = self.instrument(uri)
        temp = self._search_json(uri, 'data_label')
        if temp is None:
            temp = self._search_fits(uri, 'DATALAB')
        return obs_file_relationship.repair_data_label(uri.split('/')[-1], temp, instrument)

    def data_size(self, uri):
        return self._search_json(uri, 'data_size')

    def dec(self, uri):
        return self._search_json(uri, 'dec')

    def detector_binning(self, uri):
        return self._search_json(uri, 'detector_binning')

    def disperser(self, uri):
        return self._search_json(uri, 'disperser')

    def exposure_time(self, uri):
        return self._search_json(uri, 'exposure_time')

    def data_md5(self, uri):
        return self._search_json(uri, 'data_md5')

    def filter_name(self, uri):
        return self._search_json(uri, 'filter_name')

    def focal_plane_mask(self, uri):
        return self._search_json(uri, 'focal_plane_mask')

    def instrument(self, uri):
        temp = self._search_fits(uri, 'INSTRUME')
        if temp is None:
            if temp is None:
                temp = self._search_json(uri, 'instrument')
        result = None
        if temp is not None:
            result = repair_instrument(temp)
        return result

    def max_exputend(self, uri):
        if self._max_exputend.get(uri) is None:
            headers = self._reader.headers.get(uri)
            if headers is not None and len(headers) > 0:
                exputend_values = deque()
                for header in headers:
                    if header.get('EXPUTEND') is not None:
                        date_obs = header.get('DATE-OBS')
                        temp = header.get('EXPUTEND')
                        exputend_values.append(f'{date_obs} {temp}')

                if len(exputend_values) >= 2:
                    start = get_datetime_mjd(mc.make_datetime(exputend_values.popleft()))
                    end = get_datetime_mjd(mc.make_datetime(exputend_values.pop()))
                    if end < start:
                        # in case the observation crosses midnight
                        end = end + TimeDelta(1.0 * units.day)
                    self._max_exputend[uri] = end.value
        return self._max_exputend.get(uri)

    def mode(self, uri):
        return self._search_json(uri, 'mode')

    def object(self, uri):
        return self._search_json(uri, 'object')

    def observation_class(self, uri):
        result = self._search_json(uri, 'observation_class')
        if result is None:
            result = self._search_fits(uri, 'OBSCLASS')
        return result

    def observation_type(self, uri):
        result = self._search_json(uri, 'observation_type')
        if result is None:
            result = self._search_fits(uri, 'OBSTYPE')
        return result

    def program_id(self, uri):
        return self._search_json(uri, 'program_id')

    def ra(self, uri):
        return self._search_json(uri, 'ra')

    def reduction(self, uri):
        return self._search_json(uri, 'reduction')

    def release(self, uri):
        return self._search_json(uri, 'release')

    def spectroscopy(self, uri):
        return self._search_json(uri, 'spectroscopy')

    def telescope(self, uri):
        return self._search_json(uri, 'telescope')

    def types(self, uri):
        return self._search_json(uri, 'types')

    def ut_datetime(self, uri):
        temp = self._search_json(uri, 'ut_datetime')
        result = None
        if temp is not None:
            result = mc.make_datetime(temp)
        return result

    def _search_json(self, uri, lookup_key):
        return self._reader.json_metadata.get(uri).get(lookup_key)

    def _search_fits(self, uri, lookup_key):
        temp = None
        headers = self._reader.headers.get(uri)
        if headers is not None and len(headers) > 0:
            # if headers are None, the file is proprietary at
            # archive.gemini.edu, and cannot be retrieved
            temp = headers[0].get(lookup_key)
            if temp is None and len(headers) > 1:
                temp = headers[1].get(lookup_key)
        return temp


class GeminiMetadataLookupStorageName(GeminiMetadataLookup):
    def __init__(self, storage_name):
        super().__init__(None)
        self._storage_name = storage_name

    def max_exputend(self, uri):
        if self._max_exputend.get(uri) is None:
            headers = self._storage_name.metadata.get(uri)
            if headers is not None and len(headers) > 0:
                exputend_values = deque()
                for header in headers:
                    if header.get('EXPUTEND') is not None:
                        date_obs = header.get('DATE-OBS')
                        temp = header.get('EXPUTEND')
                        exputend_values.append(f'{date_obs} {temp}')

                if len(exputend_values) >= 2:
                    start = get_datetime_mjd(mc.make_datetime(exputend_values.popleft()))
                    end = get_datetime_mjd(mc.make_datetime(exputend_values.pop()))
                    if end < start:
                        # in case the observation crosses midnight
                        end = end + TimeDelta(1.0 * units.day)
                    self._max_exputend[uri] = end.value
        return self._max_exputend.get(uri)

    def _search_json(self, uri, lookup_key):
        return self._storage_name.json_metadata.get(uri).get(lookup_key)

    def _search_fits(self, uri, lookup_key):
        temp = None
        headers = self._storage_name.metadata.get(uri)
        if headers is not None and len(headers) > 0:
            # if headers are None, the file is proprietary at
            # archive.gemini.edu, and cannot be retrieved
            temp = headers[0].get(lookup_key)
            if temp is None and len(headers) > 1:
                temp = headers[1].get(lookup_key)
        return temp


class ProvenanceFinder:
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

    def __init__(self, clients, config):
        self._use_local_files = config.use_local_files
        self._connected = mc.TaskType.SCRAPE not in config.task_types
        self._tap_client = clients.query_client
        self._gemini_session = clients.gemini_session
        self._data_sources = config.data_sources
        self._logger = logging.getLogger(self.__class__.__name__)

    def _check_caom2(self, uri, collection):
        self._logger.debug(f'Begin _check_caom2 for {uri}')
        # The URI comparison is a LIKE, because when searching for provenance,
        # sometimes there's only file ids in the file headers, and guessing
        # the extension for gemini is not predictive.
        query_string = f"""
        SELECT O.observationID, O.instrument_name
        FROM caom2.Observation AS O
        JOIN caom2.Plane AS P on P.obsID = O.obsID
        JOIN caom2.Artifact AS A on A.planeID = P.planeID
        WHERE A.uri LIKE '{uri}%'
        AND O.collection = '{collection}'
        """
        table = clc.query_tap_client(query_string, self._tap_client)
        result = None
        instrument = None
        if len(table) == 1:
            result = table[0]['observationID']
            instrument = table[0]['instrument_name']
            self._logger.debug(f'Found observation ID {result} for {uri} at CADC.')
        self._logger.debug('End _check_caom2')
        return result, instrument

    def _check_local(self, f_name):
        self._logger.debug(f'Begin _check_local for {f_name}')
        file_id = obs_file_relationship.remove_extensions(f_name)
        try_these = [
            f'{file_id}.fits',
            f'{file_id}.fits.header',
            f'{file_id}.fits.bz2',
            f'{file_id}.fits.gz',
        ]
        result = None
        instrument = None
        for data_source in self._data_sources:
            for f_name in try_these:
                fqn = path.join(data_source, f_name)
                if path.exists(fqn):
                    headers = data_util.get_local_file_headers(fqn)
                    temp = headers[0].get('DATALAB').upper()
                    if temp is not None:
                        result = headers[0].get('DATALAB')
                        instrument = headers[0].get('INSTRUME')
                        self._logger.debug(f'Found observation ID {result} for {f_name} on ' f'disk.')
                        break
        self._logger.debug('End _check_local')
        return result, instrument

    def _check_remote(self, uri):
        self._logger.debug(f'Begin _check_remote for {uri}')
        result = None
        json = retrieve_json(uri, self._logger, self._gemini_session)
        f_id = obs_file_relationship.remove_extensions(uri.split('/')[-1])
        for ii in json:
            y = obs_file_relationship.remove_extensions(ii.get('name'))
            if y == f_id:
                result = ii.get('data_label')
                instrument = ii.get('instrument')
                break
        self._logger.debug(f'End _check_remote with result {result}')
        return result, instrument

    def get(self, uri):
        """
        :param uri: Artifact URI at CADC
        """
        ignore_scheme, collection, f_name = mc.decompose_uri(uri)
        if self._connected:
            result = None
            instrument = None
            if self._use_local_files:
                result, instrument = self._check_local(f_name)
            if result is None:
                result, instrument = self._check_caom2(uri, collection)
            if result is None:
                result, instrument = self._check_remote(uri)
        else:
            result, instrument = self._check_local(f_name)
        repaired_data_label = None
        if result is not None:
            repaired_data_label = obs_file_relationship.repair_data_label(f_name, result, instrument)
        return repaired_data_label


class GeminiMetaVisitRunnerMeta(MetaVisitRunnerMeta):
    """
    Defines the pipeline step for Collection creation or augmentation by a visitor of metadata into CAOM.
    """

    def __init__(
        self,
        clients,
        config,
        meta_visitors,
        reporter,
    ):
        super().__init__(
            clients=clients,
            config=config,
            meta_visitors=meta_visitors,
            reporter=reporter,
        )

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_uri}')
        # ask archive.gemini.edu for the information
        for index, source_name in enumerate(self._storage_name.source_names):
            uri = self._storage_name.destination_uris[index]
            if uri not in self._storage_name.metadata:
                self._storage_name.metadata[uri] = []
                if '.fits' in uri:
                    self._storage_name._metadata[uri] = retrieve_headers(
                        source_name, self._logger, self._clients, self._config
                    )
            # TODO - is there a time when not needing archive.gemini.edu is possible?
            if uri not in self._storage_name.json_metadata:
                json_record = retrieve_json(source_name, self._logger, self._clients.gemini_session)
                # json is an array of dicts, one dict per file, find the right dict
                for jj in json_record:
                    # choose this key, and the comparison, because the lhs can be
                    # a file id
                    f_name = uri.split('/')[-1]
                    if f_name in jj.get('filename'):
                        self._storage_name._json_metadata[uri] = jj
                        self._logger.debug(f'Adding JSON metadata for {uri}')
                        break

                if uri not in self._storage_name.file_info:
                    self._storage_name._file_info[uri] = FileInfo(
                        id=uri,
                        size=self._storage_name._json_metadata[uri].get('data_size'),
                        name=self._storage_name._json_metadata[uri].get('filename'),
                        md5sum=self._storage_name._json_metadata[uri].get('data_md5'),
                        lastmod=mc.make_datetime(self._storage_name._json_metadata[uri].get('lastmod')),
                        file_type=data_util.get_file_type(self._storage_name._json_metadata[uri].get('filename')),
                        encoding=data_util.get_file_encoding(self._storage_name._json_metadata[uri].get('filename')),
                    )
            if uri not in self._storage_name.file_info:
                if self._config.use_local_files:
                    self._storage_name._file_info[uri] = data_util.get_local_file_info(source_name)
                else:
                    # archive.gemini.edu lookup failed
                    self._storage_name._file_info[uri] = self._clients.data_client.info(uri)

        if not self._storage_name.obs_id:
            self._storage_name.obs_id = obs_file_relationship.repair_data_label(
                self._storage_name.file_name,
                self._storage_name._json_metadata.get(self._storage_name.file_uri).get('data_label'),
                self._storage_name._json_metadata.get(self._storage_name.file_uri).get('instrument'),
            )

        self._logger.debug('End _set_preconditions')


class GeminiOrganizeExecutesRunnerMeta(OrganizeExecutesRunnerMeta):
    """A class that extends OrganizeExecutes to handle the choosing of the correct executors based on the config.yml.
    Attributes:
        _needs_delete (bool): if True, the CAOM repo action is delete/create instead of update.
        _reporter: An instance responsible for reporting the execution status.
    Methods:
        _choose():
            Determines which descendants of CaomExecute to instantiate based on the content of the config.yml
            file for an application.
    """

    def __init__(
        self,
        config,
        meta_visitors,
        data_visitors,
        needs_delete=False,
        store_transfer=None,
        modify_transfer=None,
        clients=None,
        reporter=None,
    ):
        super().__init__(
            config,
            meta_visitors,
            data_visitors,
            store_transfer=store_transfer,
            modify_transfer=modify_transfer,
            clients=clients,
            reporter=reporter,
            needs_delete=needs_delete,
        )

    def _choose(self):
        """The logic that decides which descendants of CaomExecute to instantiate. This is based on the content of
        the config.yml file for an application.
        """
        super()._choose()
        for task_type in self.task_types:
            if task_type == mc.TaskType.INGEST or task_type == mc.TaskType.VISIT:
                if self._needs_delete:
                    raise mc.CadcException('No need identified for this yet.')
                else:
                    self._logger.debug(f'Choosing executor GeminiMetaVisit for {task_type}.')
                    self._executors = []  # over-ride the default choice.
                    self._executors.append(
                        GeminiMetaVisitRunnerMeta(self._clients, self.config, self._meta_visitors, self._reporter)
                    )

    def can_use_single_visit(self):
        return False


def repair_instrument(in_name):
    if in_name == 'ALOPEKE':
        # because the value in JSON is a different case than the value in
        # the FITS header
        in_name = 'Alopeke'
    elif in_name == 'ZORRO':
        in_name = 'Zorro'
    elif in_name == 'oscir':
        # for these observations:
        # GN-2001A-C-16-3-016
        # GN-2001A-C-2-14-015
        # GN-2001A-C-2-2-002
        # GN-2001A-C-2-3-003
        # GN-2001A-C-2-4-004
        # GN-2001A-C-2-5-005
        # GN-2001A-C-2-6-006
        # GN-2001A-C-2-7-007
        # GN-2001A-C-2-8-009
        # GN-2001A-C-2-9-010
        in_name = 'OSCIR'
    return Inst(in_name)


def retrieve_headers(source_name, logger, clients, config):
    result = None
    if config.use_local_files:
        result = data_util.get_local_file_headers(source_name)
    else:
        try:
            result = clients.data_client.get_head(f'{config.scheme}:{config.collection}/{path.basename(source_name)}')
        except exceptions.UnexpectedException as _:
            # the exceptions.NotFoundException is turned into exceptions.UnexpectedException in data_util
            # the header is not at CADC, retrieve it from archive.gemini.edu
            result = retrieve_gemini_headers(source_name, logger, clients.gemini_session)
    return result


def retrieve_gemini_headers(source_name, logger, session):
    logger.debug(f'Begin retrieve_gemini_headers for {source_name}')
    header_url = f'{HEADER_URL}{source_name}'
    # Open the URL and fetch the JSON document for the observation
    response = None
    try:
        response = session.get(header_url, timeout=20)
        response.raise_for_status()
        headers = data_util.make_headers_from_string(response.text)
    finally:
        if response is not None:
            response.close()
    logger.debug(f'End retrieve_headers')
    return headers


def retrieve_json(source_name, logger, session):
    # source name is a file id, because that's the only part that's
    # required to be unique for a retrieval from archive.gemini.edu
    #
    logger.debug(f'Begin retrieve_json for {source_name}')
    gemini_url = f'{GEMINI_METADATA_URL}{source_name}'
    # Open the URL and fetch the JSON document for the observation
    response = None
    try:
        response = mc.query_endpoint_session(gemini_url, session)
        metadata = response.json()
    finally:
        if response is not None:
            response.close()
    if len(metadata) == 0:
        raise mc.CadcException(f'Could not find JSON record for {source_name} at ' f'{gemini_url}.')
    logger.debug(f'End retrieve_json')
    return metadata
