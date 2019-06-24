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

from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc

import gem2caom2.external_metadata as em
from gem2caom2 import gem_name

__all__ = ['TapNoPreviewQuery', 'ObsFileRelationshipQuery']


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

        :param config ManageComposable.Config for gem2caom2
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
                "  AND A.uri like '{}:{}%fits' " \
                "  GROUP BY A.planeID " \
                "  HAVING COUNT(A.artifactID) = 1 ) " \
                "AND O.maxLastModified >= '{}' " \
                "AND O.maxLastModified < '{}' " \
                "AND P.dataRelease <= '{}' " \
                "ORDER BY O.maxLastModified ASC " \
                "".format(self.config.collection, gem_name.SCHEME,
                          self.config.archive, prev_exec_date, exec_date,
                          self.max_ts)
        result = ac.query_tap(query, self.config)
        return [ii for ii in result['observationID']]


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
        return em.get_gofr().subset(prev_exec_date, exec_date)
