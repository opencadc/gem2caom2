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

from datetime import datetime

from astropy.table import Table

from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc

import gem2caom2.external_metadata as em
import gem2caom2.gem_name as gem_name
import gem2caom2.obs_file_relationship as ofr

__all__ = ['TapNoPreviewQuery', 'ObsFileRelationshipQuery',
           'ArchiveGeminiEduQuery', 'EduQueryFilePre']

# See the definition of 'canonical' here, for why it matters in the URL:
# https://archive.gemini.edu/help/api.html
# ssummary = smallest query result containing file names and data labels
GEMINI_SSUMMARY_DATA = \
    'https://archive.gemini.edu/ssummary/notengineering/notFail/canonical'


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
        max_ts_s = em.get_gofr().get_max_timestamp()
        super(ObsFileRelationshipQuery, self).__init__(max_ts_s)

    def todo(self, prev_exec_date, exec_date):
        """
        Get the set of observation IDs from the in-memory storage of the
        file from Paul. The results are chunked by timestamps.

        :param prev_exec_date datetime start of the timestamp chunk
        :param exec_date datetime end of the timestamp chunk
        :return: a list of CAOM Observation IDs.
        """
        subset = em.get_gofr().subset(prev_exec_date, exec_date)
        # subset entries look like:
        # GEMINI OBSID TIMESTAMP
        # so extract the OBSID value
        temp = [ii.split()[1] for ii in subset]
        obs_ids = list(set(temp))
        return obs_ids

    def initialize(self):
        """Do nothing."""
        pass


class ArchiveGeminiEduQuery(mc.Work):

    def __init__(self, max_ts):
        super(ArchiveGeminiEduQuery, self).__init__(max_ts.timestamp())
        self._current_work_list = None
        self._encountered_max_records = False

    def todo(self, prev_exec_date, exec_date):
        """
        Get the set of observation IDs by querying archive.gemini.edu. The
        results are chunked by timestamps.

        :param prev_exec_date datetime start of the timestamp chunk
        :param exec_date datetime end of the timestamp chunk
        :return: a list of GemName instances, where GemName is an
            extension of ec.StorageName
        """
        # queries that includes times are not supported, so strip that
        # information from the query URL
        date_str = datetime.strftime(prev_exec_date, '%Y%m%d')
        obs_ids = {}
        ssummary_url = f'{GEMINI_SSUMMARY_DATA}/{date_str}/'
        logging.debug('Querying {}'.format(ssummary_url))
        response = None
        try:
            response = mc.query_endpoint(ssummary_url)
            if response is None:
                logging.warning(
                    'Could not query {}'.format(ssummary_url))
            else:
                obs_ids, max_date = self.parse_ssummary_page(
                    response.text, prev_exec_date, exec_date)
                response.close()
        finally:
            if response is not None:
                response.close()

        self._current_work_list = obs_ids
        # unique-ify the GemName instances
        return [obs_ids[ii] for ii in obs_ids]

    def initialize(self):
        pofr = ofr.PartialObsFileRelationship(
            self._current_work_list, self.max_ts_s)
        em.set_ofr(pofr)

    def check_max_records(self):
        if self._encountered_max_records:
            # retrieved the maximum number of rows - need to try to back-fill
            raise mc.CadcException(
                'Retrieved the maximum number of query rows from '
                'archive.gemini.edu. Run gem_run_edu_filepre_query.')

    def parse_ssummary_page(self, html_string, start_date, end_date):
        """Parse the html returned from archive.gemini.edu.

        :param html_string - response.text
        :param start_date datetime.datetime
        :param end_date datetime.datetime
        :return list of GemName instances for processing, plus the max
            last modified time of the files associated with the StorageName
            instances. A GemName instance captures the relationship
            between a file name and an observation ID.
        """
        work_list = {}
        max_date = start_date
        # column 0 == file name
        # column 1 == data label
        # column 2 == last modified
        if (html_string is not None and len(html_string) > 0 and
                'Your search returned no results' not in html_string):
            temp = Table.read(html_string, format='html')
            work_list_array = temp.as_array(
                names=[temp.colnames[0], temp.colnames[1], temp.colnames[2]])

            # make StorageName instances
            for ii in work_list_array:
                file_id = gem_name.GemName.remove_extensions(
                    ii[0].replace('-md!', '').replace('-fits!', ''))
                file_time = datetime.strptime(ii[2], '%Y-%m-%d %H:%M:%S')
                if start_date < file_time <= end_date:
                    logging.debug(
                        f'Adding {file_id} with timestamp {file_time} to list.')
                    storage_name = gem_name.GemName(
                        obs_id=ii[1], file_name=f'{file_id}.fits',
                        file_id=file_id)
                    work_list[file_id] = storage_name
                    max_date = max(max_date, file_time)

        if (html_string is not None and
                'search generated more than the limit of 2500' in html_string):
            self._encountered_max_records = True
        return work_list, max_date


class EduQueryFilePre(mc.Work):

    def __init__(self, max_ts):
        super(EduQueryFilePre, self).__init__(max_ts.timestamp())
        self._current_work_list = None

    def todo(self, prev_exec_date, exec_date):
        """
        Get the set of observation IDs by querying archive.gemini.edu. The
        results are chunked by timestamps. This class also chunks by
        filename prefix, in case the original query for the day exceeds
        the current limit of records retrieved (2500).

        :param prev_exec_date datetime start of the timestamp chunk
        :param exec_date datetime end of the timestamp chunk
        :return: a list of GemName instances, where GemName is an
            extension of ec.StorageName
        """
        # queries that includes times are not supported, so strip that
        # information from the query URL
        date_str = datetime.strftime(prev_exec_date, '%Y%m%d')
        obs_ids = {}
        for ii in range(0, 10):
            for prefix in ['S', 'N']:
                ssummary_url = f'{GEMINI_SSUMMARY_DATA}/{date_str}/' \
                               f'filepre={prefix}{date_str}S{ii}'
                logging.debug('Querying {}'.format(ssummary_url))
                response = None
                try:
                    response = mc.query_endpoint(ssummary_url)
                    if response is None:
                        logging.warning(
                            'Could not query {}'.format(ssummary_url))
                    else:
                        temp, max_date = self.parse_ssummary_page(
                                response.text, prev_exec_date, exec_date)
                        response.close()
                        temp_obs_ids = obs_ids
                        obs_ids = {**temp, **temp_obs_ids}
                finally:
                    if response is not None:
                        response.close()

        self._current_work_list = obs_ids
        # unique-ify the GemName instances
        temp = list(set(obs_ids.keys()))
        return [obs_ids[ii] for ii in temp]
