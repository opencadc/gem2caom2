# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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

from collections import OrderedDict
from datetime import datetime

from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc

from gem2caom2 import scrape, external_metadata

__all__ = ['TapNoPreviewQuery', 'TapRecentlyPublicQuery',
           'ObsFileRelationshipQuery', 'FileListingQuery']

# See the definition of 'canonical' here, for why it matters in the URL:
# https://archive.gemini.edu/help/api.html
# ssummary = smallest query result containing file names and data labels
GEMINI_SSUMMARY_DATA = \
    'https://archive.gemini.edu/ssummary/notengineering/NotFail/canonical'


class TapNoPreviewQuery(mc.Work):

    def __init__(self, max_ts, config):
        # using max_ts_s as the proprietary date - time when the data
        # goes public
        #
        # limit the query for this, since the point is to work with previews
        # that are only available after the proprietary period is complete
        #
        super(TapNoPreviewQuery, self).__init__(max_ts.timestamp())
        self.config = config
        self.max_ts = max_ts  # type datetime

    def todo(self, prev_exec_date, exec_date):
        """
        Get the set of observation IDs that do not have preview or
        thumbnail artifacts. The results are chunked by timestamps, and
        limited to those entries that are public.

        :param prev_exec_date datetime start of the timestamp chunk
        :param exec_date datetime end of the timestamp chunk
        :return: a list of CAOM Observation IDs.
        """
        logging.debug('Entering todo')
        query = "SELECT O.observationID " \
                "FROM caom2.Observation AS O " \
                "JOIN caom2.Plane AS P ON O.obsID = P.obsID " \
                "WHERE P.planeID IN ( " \
                "  SELECT A.planeID " \
                "  FROM caom2.Observation AS O " \
                "  JOIN caom2.Plane AS P ON O.obsID = P.obsID " \
                "  JOIN caom2.Artifact AS A ON P.planeID = A.planeID " \
                "  WHERE O.collection = '{}' " \
                "  GROUP BY A.planeID " \
                "  HAVING COUNT(A.artifactID) = 1 ) " \
                "AND O.maxLastModified >= '{}' " \
                "AND O.maxLastModified < '{}' " \
                "AND P.dataRelease <= '{}' " \
                "ORDER BY O.maxLastModified ASC " \
                "".format(self.config.collection, prev_exec_date, exec_date,
                          self.max_ts)
        result = ac.query_tap(query, self.config)
        return [ii for ii in result['observationID']]

    def initialize(self):
        """Do nothing."""
        pass


class TapRecentlyPublicQuery(mc.Work):

    def __init__(self, max_ts, config):
        super(TapRecentlyPublicQuery, self).__init__(max_ts.timestamp())
        self.config = config
        self.max_ts = max_ts  # type datetime

    def todo(self, prev_exec_date, exec_date):
        """
        Get the set of observation IDs that do not have preview or
        thumbnail artifacts. Limit the entries by time-boxing on
        dataRelease.

        :param prev_exec_date datetime start of the timestamp chunk
        :param exec_date datetime end of the timestamp chunk
        :return: a list of CAOM Observation IDs.
        """

        logging.debug('Entering todo')
        query = "SELECT O.observationID " \
                "FROM caom2.Observation AS O " \
                "JOIN caom2.Plane AS P ON O.obsID = P.obsID " \
                "WHERE P.planeID IN ( " \
                "  SELECT A.planeID " \
                "  FROM caom2.Observation AS O " \
                "  JOIN caom2.Plane AS P ON O.obsID = P.obsID " \
                "  JOIN caom2.Artifact AS A ON P.planeID = A.planeID " \
                "  WHERE O.collection = '{}' " \
                "  GROUP BY A.planeID " \
                "  HAVING COUNT(A.artifactID) = 1 ) " \
                "AND P.dataRelease > '{}' " \
                "AND P.dataRelease <= '{}' " \
                "ORDER BY O.maxLastModified ASC " \
                "".format(self.config.collection, prev_exec_date, exec_date)
        result = ac.query_tap(query, self.config)
        return [ii for ii in result['observationID']]

    def initialize(self):
        """Do nothing."""
        pass


class ObsFileRelationshipQuery(mc.Work):

    def __init__(self):
        max_ts_s = external_metadata.get_gofr().get_max_timestamp()
        super(ObsFileRelationshipQuery, self).__init__(max_ts_s)

    def todo(self, prev_exec_date, exec_date):
        """
        Get the set of observation IDs from the in-memory storage of the
        file from Paul. The results are chunked by timestamps.

        :param prev_exec_date datetime start of the timestamp chunk
        :param exec_date datetime end of the timestamp chunk
        :return: a list of CAOM Observation IDs.
        """
        subset = external_metadata.get_gofr().subset(prev_exec_date, exec_date)
        # subset entries look like:
        # GEMINI OBSID TIMESTAMP
        # so extract the OBSID value
        temp = [ii.split()[1] for ii in subset]
        obs_ids = list(set(temp))
        return obs_ids

    def initialize(self):
        """Do nothing."""
        pass


class FileListingQuery(mc.Work):
    """
    Get the set of file ids to process by querying archive.gemini.edu.

    Time-based querying of that site is limited to either a single UT Date, or
    a UT Date Range.

    The site also modifies files during post-processing, so querying a
    particular UT Date at one point in time is not guaranteed to return the
    same results as querying that same UT Date at a later point in time.

    This class attempts to deal with these two issues by querying from
    'now - 14 days' to 'now', and identifying for processing the files that
    have last modified timestamps greater than the bookmarked last processed
    time, as stored in the state file.
    """

    def __init__(self, last_processed_time):
        # last_processed_time is a datetime.datetime
        self._last_processed_time = last_processed_time
        self._max_records_encountered = False
        self._max_ts_s = datetime.utcnow().timestamp()

    def todo(self):
        # now decide on a timestamp for when to start ingestion
        # for this attempt, say 14 days prior to the last ingestion
        # timestamp - WAG
        now = datetime.utcnow()
        start_time_ts = self._last_processed_time.timestamp() - 14 * 1440 * 60

        # this list of work is the query for the day, where there are no time
        # considerations, because they is no accurate timestamp information
        # to be had from the summary pages

        prev_exec_time_ts = start_time_ts
        # work in intervals of days, because those are the units of query from
        # archive.gemini.edu
        exec_time_ts = min(
            mc.increment_time(prev_exec_time_ts, 1440).timestamp(),
            self._max_ts_s)

        entries = {}  # dict - keys are last modified, file names are value
        logging.info(f'Querying archive.gemini.edu from '
                     f'{self._last_processed_time} to {now}')
        if prev_exec_time_ts < self._max_ts_s:
            while exec_time_ts <= self._max_ts_s:
                logging.debug(f'Checking time-box for day '
                              f'{datetime.fromtimestamp(prev_exec_time_ts)}')
                temp = scrape.read_json_file_list_page(
                    prev_exec_time_ts, self._last_processed_time.timestamp())
                if len(temp) == 2500:
                    self._max_records_encountered = True
                temp_entries = entries
                entries = {**temp, **temp_entries}

                if exec_time_ts == self._max_ts_s:
                    break

                prev_exec_time_ts = exec_time_ts
                exec_time_ts = min(
                    mc.increment_time(prev_exec_time_ts, 1440).timestamp(),
                    self._max_ts_s)

                # TODO - size and metadata checks, to see if those have
                # changed, before adding an entry to the todo list
                # no work has been done yet, so don't save any state anywhere

        ordered_file_names = OrderedDict(sorted(entries.items()))
        logging.info(f'Found {len(ordered_file_names)} entries to process.')
        return ordered_file_names
