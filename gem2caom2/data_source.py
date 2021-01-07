# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
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

from collections import deque
from datetime import datetime, timezone

from caom2pipe import data_source_composable as dsc
from caom2pipe import manage_composable as mc
from gem2caom2.scrape import read_json_file_list_page


__all__ = ['FileListIncrementalSource', 'GEM_BOOKMARK', 'IncrementalSource',
           'PublicIncremental']

GEM_BOOKMARK = 'gemini_timestamp'


class IncrementalSource(dsc.DataSource):
    """Implements the identification of the work to be done, by querying
    archive.gemini.edu's incremental endpoint, in time-boxed chunks."""

    def __init__(self):
        super(IncrementalSource, self).__init__(config=None)
        self._logger = logging.getLogger(__name__)

    def get_time_box_work(self, prev_exec_time, exec_time):
        """
        :param prev_exec_time datetime start of the timestamp chunk
        :param exec_time datetime end of the timestamp chunk
        :return: a list of file names with time they were modified from
            archive.gemini.edu, structured as an astropy Table (for now).
        """

        self._logger.debug('Entering todo')
        # datetime format 2019-12-01T00:00:00.000000
        prev_dt_str = prev_exec_time.strftime(mc.ISO_8601_FORMAT)
        exec_dt_str = exec_time.strftime(mc.ISO_8601_FORMAT)
        url = f'http://arcdev.gemini.edu/jsonsummary/sr=179/notengineering/' \
              f'NotFail/entrytime={prev_dt_str}%20{exec_dt_str}?order_by=' \
              f'entrytime'

        # needs to be ordered by timestamps when processed
        entries = deque()
        response = None
        try:
            response = mc.query_endpoint(url)
            if response is None:
                logging.warning(f'Could not query {url}.')
            else:
                metadata = response.json()
                response.close()
                if metadata is not None:
                    if len(metadata) == 0:
                        self._logger.warning(f'No query results returned for '
                                             f'interval from {prev_exec_time} '
                                             f'to {exec_time}.')
                    else:
                        for entry in metadata:
                            file_name = entry.get('name')
                            entries.append(dsc.StateRunnerMeta(
                                file_name, entry.get('lastmod')))
        finally:
            if response is not None:
                response.close()
        return entries


class FileListIncrementalSource(dsc.DataSource):
    """Implements the identification of the work to be done, by querying
    archive.gemini.edu's jsonfilelist endpoint, in time-boxed chunks."""

    def __init__(self, config):
        super(FileListIncrementalSource, self).__init__(config=None)
        self._logger = logging.getLogger(__name__)
        state = mc.State(config.state_fqn)
        temp = state.get_bookmark(GEM_BOOKMARK)
        # make sure last_processed_time is type float
        self._last_processed_time = mc.increment_time(temp, 0).timestamp()
        # now decide on a timestamp for when to start ingestion
        # for this attempt, say 14 days prior to the last ingestion
        # timestamp - WAG
        self._start_time_ts = self._last_processed_time - 14 * 1440 * 60
        self._max_records_encountered = False

    def get_time_box_work(self, prev_exec_time, exec_time):
        """
        :param prev_exec_time datetime start of the timestamp chunk
        :param exec_time datetime end of the timestamp chunk
        :return: a list of file names with time they were modified from
            archive.gemini.edu, structured as an astropy Table (for now).
        """
        self._logger.debug(f'Begin get_time_box_work from '
                           f'{prev_exec_time} to {exec_time}.')
        temp = read_json_file_list_page(
            prev_exec_time, self._last_processed_time)
        if len(temp) == 2500:
            self._max_records_encountered = True

        entries = deque()
        for key, value in temp.items():
            entries.append(dsc.StateRunnerMeta(
                value, datetime.fromtimestamp(key).isoformat()))
        self._logger.debug('End get_time_box_work.')
        return entries

    @property
    def max_records_encountered(self):
        return self._max_records_encountered


class PublicIncremental(dsc.QueryTimeBoxDataSource):
    """Implements the identification of the work to be done, by querying
    the local TAP service for files that have recently gone public."""

    def __init__(self, config):
        super(PublicIncremental, self).__init__(config)
        self._logger = logging.getLogger(__name__)

    def get_time_box_work(self, prev_exec_time, exec_time):
        """
        :param prev_exec_time datetime start of the timestamp chunk
        :param exec_time datetime end of the timestamp chunk
        :return: a list of file names with time they were modified in /ams,
            structured as an astropy Table (for now).
        """

        self._logger.debug('Entering get_time_box_work')
        # datetime format 2019-12-01T00:00:00.000000
        prev_dt_str = datetime.fromtimestamp(
            prev_exec_time, tz=timezone.utc).strftime(mc.ISO_8601_FORMAT)
        exec_dt_str = datetime.fromtimestamp(
            exec_time, tz=timezone.utc).strftime(mc.ISO_8601_FORMAT)
        query = f"SELECT A.uri, A.lastModified " \
                f"FROM caom2.Observation AS O " \
                f"JOIN caom2.Plane AS P ON O.obsID = P.obsID " \
                f"JOIN caom2.Artifact AS A ON P.planeID = A.planeID " \
                f"WHERE P.planeID IN ( " \
                f"  SELECT A.planeID " \
                f"  FROM caom2.Observation AS O " \
                f"  JOIN caom2.Plane AS P ON O.obsID = P.obsID " \
                f"  JOIN caom2.Artifact AS A ON P.planeID = A.planeID " \
                f"  WHERE O.collection = '{self._config.collection}' " \
                f"  GROUP BY A.planeID " \
                f"  HAVING COUNT(A.artifactID) = 1 ) " \
                f"AND P.dataRelease > '{prev_dt_str}' " \
                f"AND P.dataRelease <= '{exec_dt_str}' " \
                f"ORDER BY O.maxLastModified ASC " \
                ""
        result = mc.query_tap_client(query, self._client)
        # results look like:
        # gemini:GEM/N20191202S0125.fits, ISO 8601

        entries = deque()
        for row in result:
            entries.append(dsc.StateRunnerMeta(
                mc.CaomName(row['uri']).file_name, row['lastModified']))
        return entries
