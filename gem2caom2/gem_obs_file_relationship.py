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
import collections
import logging

from datetime import datetime
from datetime import timezone

from caom2pipe import manage_composable as mc


class GemObsFileRelationship(object):
    """A class to hold and access the content of the observation ID/file id
    information that is provided to CADC from Gemini.

    It's made into a class, because the information is useful from both
    the gem2caom2 repo, for identifying provenance relationships, and from
    the gemHarvester2Caom2, for supporting list_observations queries.
    """

    def __init__(self, file_name):
        self.id_list = {}
        self.time_list = {}
        self.name_list = {}
        self.logger = logging.getLogger('GemObsFileRelationship')
        self._initialize_content(file_name)

    def _initialize_content(self, fqn):
        """Initialize the internal data structures that represents the
        query list from the Gemini Science Archive.

        observation_list structure: a dict, keys are last modified time,
            values are a set of observation IDs with that last modified time

        observation_id_list structure: a dict, keys are observation ID,
            values are a set of associated file names
        """
        result = self._read_file(fqn)
        # result row structure:
        # 0 = data label
        # 1 = timestamp
        # 2 = file name
        temp_content = {}
        for ii in result:
            # re-organize to be able to answer list_observations queries
            ol_key = self._make_seconds(ii[1])
            if ol_key in temp_content:
                if ii[0] not in temp_content[ol_key]:
                    temp_content[ol_key].append(ii[0])
            else:
                temp_content[ol_key] = [ii[0]]
            # re-organize to be able to answer get_observation queries
            if ii[0] in self.id_list:
                self.id_list[ii[0]].append(ii[2])
            else:
                self.id_list[ii[0]] = [ii[2]]
            if ii[2] in self.name_list:
                self.name_list[ii[2]].append(ii[0])
            else:
                self.name_list[ii[2]] = [ii[0]]

        # this structure means an observation ID occurs more than once with
        # different last modified times
        self.time_list = collections.OrderedDict(sorted(temp_content.items(),
                                                        key=lambda t: t[0]))
        self.logger.debug('Observation list initialized in memory.')

    def _read_file(self, fqn):
        """Read the .txt file from Gemini, and make it prettier ...
        where prettier means stripping whitespace, query output text, and
        making an ISO 8601 timestamp from something that looks like this:
        ' 2018-12-17 18:19:27.334144+00 '

        or this:
        ' 2018-12-17 18:19:27+00 '

        :return a list of lists, where the inner list consists of an
            observation ID, a last modified date/time, and a file name.

        File structure indexes:
        0 == data label
        1 == file name
        3 == last modified date/time
        """
        results = []
        try:
            with open(fqn) as f:
                for row in f:
                    temp = row.split('|')
                    if len(temp) > 1 and 'data_label' not in row:
                        time_string = temp[3].strip().replace(' ', 'T')
                        if len(temp[0].strip()) > 1:
                            results.append(
                                [temp[0].strip(), time_string, temp[1].strip()])
                        else:
                            # no data label in the file, so use the file name
                            results.append(
                                [temp[1].strip(), time_string, temp[1].strip()])

        except Exception as e:
            self.logger.error('Could not read from csv file {}'.format(fqn))
            raise mc.CadcException(e)
        return results

    @staticmethod
    def _make_seconds(from_time):
        """Deal with the different time formats in the Gemini-supplied file
        to get the number of seconds since the epoch, to serve as an
        ordering index for the list of observation IDs.

        The obs id file has the timezone information as +00, strip that for
        returned results.
        """
        index = from_time.index('+00')
        try:
            seconds_since_epoch = datetime.strptime(from_time[:index],
                                                    '%Y-%m-%dT%H:%M:%S.%f')
        except ValueError as e:
            seconds_since_epoch = datetime.strptime(from_time[:index],
                                                    '%Y-%m-%dT%H:%M:%S')
        return seconds_since_epoch.timestamp()

    def subset(self, start=None, end=None, maxrec=None):
        if start is not None and end is not None:
            temp = self._subset(start.timestamp(), end.timestamp())
        elif start is not None:
            temp = self._subset(start.timestamp(), datetime.now().timestamp())
        elif end is not None:
            temp = self._subset(0, end.timestamp())
        else:
            temp = self._subset(0, datetime.now().timestamp())
        if maxrec is not None:
            temp = temp[:maxrec]
        return temp

    def _subset(self, start_s, end_s):
        """Get only part of the observation list, limited by timestamps."""
        self.logger.debug('Timestamp endpoints are between {} and {}.'.format(
            start_s, end_s))
        temp = []
        for ii in self.time_list:
            if start_s <= ii <= end_s:
                for jj in self.time_list[ii]:
                    temp.append(
                        '{},{}'.format(jj, datetime.fromtimestamp(ii, timezone.utc).isoformat()))
            if ii > end_s:
                break
        return temp

    def get_file_names(self, obs_id):
        if obs_id in self.id_list:
            return self.id_list[obs_id]
        else:
            return None

    def get_obs_id(self, file_name):
        if file_name in self.name_list:
            return self.name_list[file_name]
        else:
            return None
