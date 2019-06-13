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
import re

from datetime import datetime
from datetime import timedelta

from caom2pipe import manage_composable as mc

from gem2caom2 import gem_name


__all__ = ['GemObsFileRelationship', 'CommandLineBits']

HEADER_URL = 'https://archive.gemini.edu/fullheader/'


class GemObsFileRelationship(object):
    """A class to hold and access the content of the observation ID/file id
    information that is provided to CADC from Gemini.

    It's made into a class, because the information is useful from both
    the gem2caom2 repo, for identifying provenance relationships, and from
    the gemHarvester2Caom2, for supporting list_observations queries.
    """

    def __init__(self, file_name):

        # id_list structure: a dict, keys are Gemini observation IDs, values
        # are a set of associated file names. This structure supports the
        # get_observation query.

        self.id_list = {}

        # time_list structure: a dict, keys are last modified time,
        # values are a set of observation IDs as specified from Gemini
        # with that last modified time. This structure supports the
        # time-bounded queries of the Harvester.

        self.time_list = {}

        # name_list structure: a dict, keys are file ids, values are Gemini
        # observation IDs. This structure supports queries by gem2caom2
        # for determining provenance information for planes and
        # observations.

        self.name_list = {}

        # structure: a dict, keys are Gemini observation IDs, values are
        # repaired observation IDs. This structure supports generation
        # of the command line for invoking gem2caom2

        self.repaired_ids = {}

        # the repaired obs id to file name lookup

        self.repaired_names = {}

        self.logger = logging.getLogger('GemObsFileRelationship')
        self._initialize_content(file_name)

    def _initialize_content(self, fqn):
        """Initialize the internal data structures that represents the
        query list from the Gemini Science Archive.
        """
        result = self._read_file(fqn)
        # result row structure:
        # 0 = data label
        # 1 = timestamp
        # 2 = file name
        temp_content = {}
        logging.error('well it started ....')
        for ii in result:
            # re-organize to be able to answer list_observations queries
            ol_key = mc.make_seconds(ii[1])
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
            file_id = gem_name.GemName.remove_extensions(ii[2])
            if file_id in self.name_list:
                self.name_list[file_id].append([ii[0], ol_key])
            else:
                self.name_list[file_id] = [[ii[0], ol_key]]

        logging.error('After the initial bits')

        self._build_repaired_lookups()

        logging.error('After the repair bits')

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
                        if '/' in temp[0]:
                            if 'MBIAS' in temp[0]:
                                temp[0] = temp[0].replace('BIAS/MBIAS/', '')
                            elif 'PETRO' in temp[0]:
                                temp[0] = temp[0].replace(
                                    '-/NET/PETROHUE/DATAFLOW/',
                                    '')
                            elif '12CD' in temp[0] or 'EXPORT/HOME' in temp[0]:
                                temp[0] = temp[0].split('/', 1)[0]
                            else:
                                logging.warning(
                                    'Mystery data label {}'.format(temp[0]))
                        elif '?' in temp[0]:
                            if 'GS-2002A-DD-1-?' in temp[0]:
                                temp[0] = temp[0].replace('?', '11')
                            else:
                                logging.warning(
                                    'Mystery data label {}'.format(temp[0]))
                        elif '"' in temp[0]:
                            temp[0] = temp[0].replace('"', '')
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
                        '{} {} {}'.format(
                            gem_name.COLLECTION, jj,
                            datetime.fromtimestamp(ii).isoformat(
                                timespec='milliseconds')))
            if ii > end_s:
                break
        return temp

    def get_file_names(self, obs_id):
        if obs_id in self.id_list:
            return self.id_list[obs_id]
        else:
            return None

    def get_obs_id(self, file_id):
        if file_id in self.name_list:
            return self.name_list[file_id][0]
        else:
            return None

    def get_timestamp(self, file_id):
        if file_id in self.name_list:
            temp = self.name_list[file_id]
            return temp[0][1]
        else:
            return timedelta()

    @staticmethod
    def is_processed(file_name):
        """Try to determine if a Gemini file is processed, based on naming
        patterns."""
        result = True
        file_id = gem_name.GemName.remove_extensions(file_name)
        if (file_id.startswith(('S', 'N', 'GN', 'GS', 'c', 'abu')) and
                file_id.endswith(
                    ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'))):
            result = False
        elif file_id.startswith(('2', '02', '01')):
            result = False
        # TEXES naming patterns
        elif file_id.startswith('TX2') and '_raw' in file_id:
            result = False
        # OSCIR file naming pattern
        elif (file_id.startswith('r') and
              re.match('r\\w{7}_\\d{3}', file_id, flags=re.ASCII) is not None):
            result = False
        return result

    def repair_data_label(self, file_id):
        """For processed files, try to provide a consistent naming pattern,
        because data labels aren't unique within Gemini, although the files
        they refer to are, and can be in different CAOM Observations.

        Take the prefixes and suffixes on the files, that indicate the type of
        processing, and append them in upper case, to the data label, for
        uniqueness.

        DB - 07-03-19
        TEXES Spectroscopy

        Some special code will be needed for datalabels/planes.  There are no
        datalabels in the FITS header.  json metadata (limited) must be
        obtained with URL like
        https://archive.gemini.edu/jsonsummary/canonical/filepre=TX20170321_flt.2507.fits.
        Use TX20170321_flt.2507 as datalabel.  But NOTE:  *raw.2507.fits and
        *red.2507.fits are two planes of the same observation. I’d suggest we
        use ‘*raw*’ as the datalabel and ‘*red*’ or ‘*raw*’ as the appropriate
        product ID’s for the science observations.  The ‘flt’ observations do
        not have a ‘red’ plane.  The json document contains ‘filename’ if
        that’s helpful at all.  The ‘red’ files do not exist for all ‘raw’
        files.
        """
        if file_id in self.name_list:
            repaired = self.name_list[file_id][0][0]
            # if the data label is missing, the file name, including
            # extensions, is treated as the data label, so get rid of .fits
            repaired = gem_name.GemName.remove_extensions(repaired)
            self.logger.debug(
                'Gemini says data label is {} for file {}'.format(repaired,
                                                                  file_id))
            if (GemObsFileRelationship.is_processed(file_id) or
                    file_id.startswith('TX2')):
                if not file_id.startswith('TX2'):
                    repaired = repaired.split('_')[0]
                prefix = GemObsFileRelationship._get_prefix(file_id)
                suffix = GemObsFileRelationship._get_suffix(file_id, repaired)
                removals = GemObsFileRelationship._get_removals(file_id,
                                                                repaired)

                if len(prefix) > 0:
                    removals = [prefix] + suffix
                else:
                    removals = removals + suffix
                self.logger.debug(
                    'repaired {} removals {}'.format(repaired, removals))
                for ii in removals:
                    # rreplace
                    temp = repaired.rsplit(ii, 1)
                    repaired = ''.join(temp)
                    temp = repaired.rsplit(ii.upper(), 1)
                    repaired = ''.join(temp)
                    repaired = repaired.rstrip('-')
                    repaired = repaired.rstrip('_')

                # DB - 18-03-19
                # Basically any ‘mfrg’, ‘mrg’, or ‘rg’ file WITHOUT ‘add’
                # in the datalabel or name is a processed version of a raw
                # file without the datalabel suffix (filename prefix)
                if 'mfrg' == prefix or 'mrg' == prefix or 'rg' == prefix:
                    if not ('add' in suffix or 'ADD' in suffix):
                        prefix = ''
                        suffix = []

                if len(prefix) > 0:
                    repaired = '{}-{}'.format(repaired, prefix.upper())
                for ii in suffix:
                    repaired = '{}-{}'.format(repaired, ii.upper())
            else:
                repaired = file_id if repaired is None else repaired

            if file_id == 'N20181217S0266':
                repaired = 'GN-2018B-Q-133-20-001'
        else:
            logging.warning(
                'File name {} not found in the Gemini list.'.format(file_id))
            repaired = file_id
        return repaired

    @staticmethod
    def _get_prefix(file_id):
        if file_id.startswith(('p', 'P')):
            if '_FLAT' in file_id or '_COMB' in file_id:
                prefix = 'P'
            else:
                prefix = ''
        elif file_id.startswith(('TX', 'ag', 'c')):
            prefix = ''
        elif -1 < file_id.find('GN') < 14:
            prefix = file_id.split('GN', 1)[0]
        elif -1 < file_id.find('N') < 14:
            prefix = file_id.split('N', 1)[0]
        elif 'GS' in file_id:
            prefix = file_id.split('GS', 1)[0]
        elif 'S' in file_id:
            prefix = file_id.split('S', 1)[0]
        else:
            logging.warning(
                'Unrecognized file_id pattern {}'.format(file_id))
            prefix = ''
        return prefix

    @staticmethod
    def _get_suffix(file_id, data_label):
        temp = []
        suffix = []
        if '-' in file_id:
            temp = file_id.split('-')[:1]
        elif '_' in file_id:
            if file_id.startswith(('p', 'P')):
                if '_FLAT' in file_id or '_COMB' in file_id:
                    temp = file_id.split('_')[2:]
            elif file_id.startswith('TX2'):
                # when the data label is the file id, fix every
                # data label except flats
                if not data_label.startswith(file_id):
                    if '_flt' in file_id.lower():
                        temp = ['flt']
            else:
                temp = file_id.split('_')[1:]
                # logging.error('get here? suffix is {} for {}'.format(suffix, file_id))
        for ii in temp:
            if re.match('[a-zA-Z]+', ii) is not None:
                suffix.append(ii)
        return suffix

    @staticmethod
    def _get_removals(file_id, repaired):
        removals = []
        if file_id.startswith('TX2'):
            # when the data label is the file id, fix every
            # data label except flats
            if repaired.startswith(file_id):
                if '_flt' not in file_id.lower():
                    removals = ['_raw', '_red', '_sum']
            else:
                if '_flt' not in file_id.lower():
                    removals = ['raw', 'red', 'sum']
        return removals

    def _build_repaired_lookups(self):
        # for each gemini observation ID, get the file names associated with
        # that observation ID
        index = 0
        for ii in self.id_list:
            file_names = self.id_list[ii]
            # for each file name, repair the obs id, add repaired obs id
            # and file name to new structure
            for file_name in file_names:
                file_id = gem_name.GemName.remove_extensions(file_name)
                repaired_obs_id = self.repair_data_label(file_id)
                temp = gem_name.GemName.remove_extensions(ii)
                self._add_repaired_element(temp, repaired_obs_id, file_id)
            index += 1

            if index % 100000 == 0:
                logging.error('got to {} {}'.format(index, ii))

        index = 0
        # for each file name, add repaired obs ids, if they're not already
        # in the list
        for file_name in self.name_list:
            for ii in self.name_list[file_name]:
                obs_id = ii[0]
                file_id = gem_name.GemName.remove_extensions(file_name)
                repaired_obs_id = self.repair_data_label(file_id)
                self._add_repaired_element(obs_id, repaired_obs_id, file_id)
            index += 1

            if index % 100000 == 0:
                logging.error('got to {} {}'.format(index, file_name))

    def _add_repaired_element(self, obs_id, repaired_obs_id, file_id):
        if obs_id in self.repaired_ids:
            if repaired_obs_id not in self.repaired_ids[obs_id]:
                self.repaired_ids[obs_id].append(repaired_obs_id)
        else:
            self.repaired_ids[obs_id] = [repaired_obs_id]
        if repaired_obs_id in self.repaired_names:
            if file_id not in self.repaired_names[repaired_obs_id]:
                self.repaired_names[repaired_obs_id].append(file_id)
        else:
            self.repaired_names[repaired_obs_id] = [file_id]

    def get_args(self, obs_id):
        if obs_id in self.repaired_ids:
            result = []
            for repaired_id in self.repaired_ids[obs_id]:
                lineage = ''
                urls = ''
                for file_id in self.repaired_names[repaired_id]:
                    # works because the file id == product id
                    lineage += mc.get_lineage(
                        gem_name.ARCHIVE, file_id, '{}.fits'.format(file_id),
                        gem_name.SCHEME)
                    urls += '{}{}.fits'.format(HEADER_URL, file_id)
                    if file_id != self.repaired_names[repaired_id][-1]:
                        lineage += ' '
                        urls += ' '
                c = CommandLineBits(
                    '{} {}'.format(gem_name.COLLECTION, repaired_id),
                    lineage, urls)
                result.append(c)
            return result
        else:
            logging.warning(
                'Could not find observation ID {} in Gemini-provided '
                'list.'.format(obs_id))
            return []


class CommandLineBits(object):
    """Convenience class to keep the bits of command-line that are
    inter-connected together."""

    def __init__(self, obs_id='', lineage='', urls=''):
        self.obs_id = obs_id
        self.lineage = lineage
        self.urls = urls

    def __str__(self):
        return '{}\n{}\n{}'.format(self.obs_id, self.lineage, self.urls)

    @property
    def obs_id(self):
        return self._obs_id

    @obs_id.setter
    def obs_id(self, value):
        self._obs_id = value

    @property
    def lineage(self):
        return self._lineage

    @lineage.setter
    def lineage(self, value):
        self._lineage = value

    @property
    def urls(self):
        return self._urls

    @urls.setter
    def urls(self, value):
        self._urls = value
