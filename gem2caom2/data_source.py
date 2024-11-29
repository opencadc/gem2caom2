# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2024.                            (c) 2024.
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

from bs4 import BeautifulSoup
from collections import defaultdict, deque, OrderedDict
from datetime import datetime

from caom2pipe import client_composable as clc
from caom2pipe import data_source_composable as dsc
from caom2pipe.manage_composable import build_uri, CaomName, ISO_8601_FORMAT, make_datetime, query_endpoint_session
from caom2pipe.manage_composable import StorageName


__all__ = ['GEM_BOOKMARK', 'IncrementalSource', 'PublicIncremental']

GEM_BOOKMARK = 'gemini_timestamp'
MAX_ENTRIES = 500

class IncrementalSource(dsc.IncrementalDataSource):
    """Implements the identification of the work to be done, by querying
    archive.gemini.edu's incremental endpoint, in time-boxed chunks.

    OO - email 21-04-21
    entrytime is when the DB record behind the JSON being displayed was
    created.
    """

    def __init__(self, config, reader):
        super().__init__(config, start_key=GEM_BOOKMARK)
        self._max_records_encountered = False
        self._encounter_start = None
        self._encounter_end = None
        self._session = reader._session
        self._metadata_reader = reader

    def _initialize_end_dt(self):
        self._end_dt = datetime.now()

    def get_time_box_work(self, prev_exec_dt, exec_dt):
        """
        :param prev_exec_dt datetime start of the time-boxed chunk
        :param exec_dt datetime end of the time-boxed chunk
        :return: a deque of file names with time their associated JSON (DB)
            records were modified from archive.gemini.edu.
        """

        self._logger.debug(f'Begin get_time_box_work from {prev_exec_dt} to {exec_dt}.')
        # datetime format 2019-12-01T00:00:00.000000
        prev_dt_str = prev_exec_dt.strftime(ISO_8601_FORMAT)
        exec_dt_str = exec_dt.strftime(ISO_8601_FORMAT)
        url = (
            f'https://archive.gemini.edu/jsonsummary/canonical/'
            f'NotFail/notengineering/'
            f'entrytimedaterange={prev_dt_str}%20{exec_dt_str}/'
            f'?orderby=entrytime'
        )

        # needs to be ordered by timestamps when processed
        self._logger.info(f'Querying {url}')
        entries = deque()
        response = None
        try:
            response = query_endpoint_session(url, self._session)
            if response is None:
                self._logger.warning(f'Could not query {url}.')
            else:
                metadata = response.json()
                response.close()
                if metadata is not None:
                    if len(metadata) == 0:
                        self._logger.warning(
                            f'No query results returned for interval from {prev_exec_dt} to {exec_dt}.'
                        )
                    else:
                        for entry in metadata:
                            file_name = entry.get('name')
                            entry_dt = make_datetime(entry.get('entrytime'))
                            entries.append(dsc.StateRunnerMeta(file_name, entry_dt))
                            uri = build_uri(StorageName.collection, file_name, StorageName.scheme)
                            # all the other cases where add_json_record is
                            # called, there's a list as input, so conform to
                            # that typing here
                            self._metadata_reader.add_json_record(uri, [entry])
                            self._metadata_reader.add_file_info_record(uri)
        finally:
            if response is not None:
                response.close()
        if len(entries) == MAX_ENTRIES:
            self._max_records_encountered = True
            self._encounter_start = prev_exec_dt
            self._encounter_end = exec_dt
        self._reporter.capture_todo(len(entries), 0, 0)
        self._logger.debug('End get_time_box_work.')
        return entries

    @property
    def max_records_encountered(self):
        if self._max_records_encountered:
            self._logger.error(
                f'Max records window {self._encounter_start} to '
                f'{self._encounter_end}.'
            )
        return self._max_records_encountered


class PublicIncremental(dsc.QueryTimeBoxDataSource):
    """Implements the identification of the work to be done, by querying
    the local TAP service for files that have recently gone public."""

    def __init__(self, config, query_client):
        super().__init__(config)
        self._query_client = query_client

    def _initialize_end_dt(self):
        self._end_dt = datetime.now()

    def get_time_box_work(self, prev_exec_dt, exec_dt):
        """
        :param prev_exec_dt datetime start of the timestamp chunk
        :param exec_dt datetime end of the timestamp chunk
        :return: a list of file names with time they were modified in /ams,
            structured as an astropy Table (for now).
        """
        self._logger.debug('Begin get_time_box_work')
        # datetime format 2019-12-01T00:00:00.000000
        prev_dt_str = prev_exec_dt.strftime(ISO_8601_FORMAT)
        exec_dt_str = exec_dt.strftime(ISO_8601_FORMAT)
        query = (
            f"SELECT A.uri, A.lastModified "
            f"FROM caom2.Observation AS O "
            f"JOIN caom2.Plane AS P ON O.obsID = P.obsID "
            f"JOIN caom2.Artifact AS A ON P.planeID = A.planeID "
            f"WHERE P.planeID IN ( "
            f"  SELECT A.planeID "
            f"  FROM caom2.Observation AS O "
            f"  JOIN caom2.Plane AS P ON O.obsID = P.obsID "
            f"  JOIN caom2.Artifact AS A ON P.planeID = A.planeID "
            f"  WHERE O.collection = '{self._config.collection}' "
            f"  GROUP BY A.planeID "
            f"  HAVING COUNT(A.artifactID) = 1 ) "
            f"AND P.dataRelease > '{prev_dt_str}' "
            f"AND P.dataRelease <= '{exec_dt_str}' "
            f"ORDER BY O.maxLastModified ASC "
            ""
        )
        result = clc.query_tap_client(query, self._query_client)
        # results look like:
        # gemini:GEM/N20191202S0125.fits, ISO 8601

        entries = deque()
        for row in result:
            entries.append(
                dsc.StateRunnerMeta(CaomName(row['uri']).file_name, make_datetime(row['lastModified']))
            )
        self._reporter.capture_todo(len(entries), 0, 0)
        self._logger.debug('End get_time_box_work')
        return entries


class IncrementalSourceDiskfiles(dsc.IncrementalDataSource):
    """Implements the identification of the work to be done, by querying archive.gemini.edu's incremental diskfiles
    endpoint, in time-boxed chunks.

    OO - email 21-04-21
    entrytime is when the DB record behind the JSON being displayed was created.
    """

    def __init__(self, config, gemini_session, storage_name_ctor, filter_cache):
        super().__init__(config, start_key=GEM_BOOKMARK)
        self._max_records_encountered = False
        self._encounter_start = None
        self._encounter_end = None
        self._session = gemini_session
        self._storage_name_ctor = storage_name_ctor
        self._filter_cache = filter_cache

    def _initialize_end_dt(self):
        self._end_dt = datetime.now()

    def get_time_box_work(self, prev_exec_dt, exec_dt):
        """
        :param prev_exec_dt datetime start of the time-boxed chunk
        :param exec_dt datetime end of the time-boxed chunk
        :return: a deque of file names with time their associated JSON (DB) records were modified from
            archive.gemini.edu.
        """

        self._logger.debug(f'Begin get_time_box_work from {prev_exec_dt} to {exec_dt}.')
        self._max_records_encountered = False
        # datetime format 2019-12-01T00:00:00 => no microseconds in the url
        prev_exec_dt_iso = prev_exec_dt.replace(microsecond=0)
        exec_dt_iso = exec_dt.replace(microsecond=0)
        url = (
            f'https://archive.gemini.edu/diskfiles/entrytimedaterange='
            f'{prev_exec_dt_iso.isoformat()}--{exec_dt_iso.isoformat()}'
        )

        # needs to be ordered by timestamps when processed
        self._logger.info(f'Querying {url}')
        entries = deque()
        response = None
        try:
            response = query_endpoint_session(url, self._session)
            if response is None:
                self._logger.warning(f'Could not query {url}.')
            else:
                # x is a dict, sorted by the key
                # key is the timestamp
                # value is the html meta that comes back, also in a dict
                metadata = self._parse_diskfiles_response(response.text)
                response.close()
                if len(metadata) == 0:
                    self._logger.warning(
                        f'No query results returned for interval from {prev_exec_dt} to {exec_dt}.'
                    )
                else:
                    for entry_dt, values  in metadata.items():
                        for value in values:
                            file_name = value.get('filename')
                            storage_name = self._storage_name_ctor(file_name, self._filter_cache)
                            storage_name.file_info[storage_name.destination_uris[0]] = value
                            entries.append(dsc.RunnerMeta(storage_name, entry_dt))
        finally:
            if response is not None:
                response.close()
        if len(entries) == MAX_ENTRIES:
            self._logger.error(f'Max records window {self._encounter_start} to {self._encounter_end}.')
            self._max_records_encountered = True
            self._encounter_start = prev_exec_dt
            self._encounter_end = exec_dt
        self._reporter.capture_todo(len(entries), 0, 0)
        self._logger.debug('End get_time_box_work.')
        return entries

    def _parse_diskfiles_response(self, html_string):
        temp = defaultdict(list)
        soup = BeautifulSoup(html_string, features='lxml')
        rows = soup.find_all('tr', {'class':'alternating'})
        for row in rows:
            cells = row.find_all('td')
            entry_time = make_datetime(cells[5].text.strip())  # what the https query is keyed on
            value = {
                'filename': cells[0].text.strip(),
                'datalabel': cells[1].text.strip(),
                'instrument': cells[2].text.strip(),
                'lastmod': make_datetime(cells[6].text.strip()),
                'data_size': cells[10].text.strip(),
                'data_md5': cells[11].text.strip(),
            }
            temp[entry_time].append(value)
        result = OrderedDict(sorted(temp.items()))
        return result

    def max_records_encountered(self):
        return self._max_records_encountered
