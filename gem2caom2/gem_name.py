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

from caom2pipe import manage_composable as mc

from gem2caom2 import external_metadata as em
from gem2caom2 import obs_file_relationship as ofr
from gem2caom2 import scrape


__all__ = ['GemName', 'COLLECTION', 'ARCHIVE', 'SCHEME', 'GemNameBuilder']


COLLECTION = 'GEMINI'
ARCHIVE = 'GEM'
SCHEME = 'gemini'


class GemName(mc.StorageName):
    """Naming rules:
    - support mixed-case file name storage, exception for extensions, and
            mixed-case obs id values - the case the inputs are provided in are
            assumed to be correct.
    - support uncompressed files in storage
    """

    GEM_NAME_PATTERN = '*'

    def __init__(self, fname_on_disk=None, file_name=None, obs_id=None,
                 file_id=None):
        logging.debug('parameters fname_on_disk {} file_name {}'
                      ' obs id {} file id {}'.format(fname_on_disk,
                                                     file_name,
                                                     obs_id,
                                                     file_id))
        # try to set the file name, if that information is available

        # file_name is assumed to be the file name in ad
        # because the GEM files are stored uncompressed,
        # while the files available from Gemini are bz2.
        self._file_name = None
        self._file_id = None
        if file_name is not None:
            self._file_id = GemName.get_file_id(file_name)
            self.file_name = file_name
        if fname_on_disk is not None:
            self._file_id = GemName.get_file_id(fname_on_disk)
            self.file_name = fname_on_disk
        if obs_id is not None:
            self._obs_id = obs_id
        if fname_on_disk is None and file_name is None and obs_id is None:
            raise mc.CadcException('Require a name.')
        super(GemName, self).__init__(
            obs_id=obs_id, collection=ARCHIVE,
            collection_pattern=GemName.GEM_NAME_PATTERN,
            fname_on_disk=self.file_name,
            scheme=SCHEME)
        if self._obs_id is None:
            temp = em.get_gofr().get_obs_id(self._file_id)
            if temp is not None:
                self._obs_id = GemName.remove_extensions(temp)
        if self._obs_id is None or self._obs_id == 'None':
            # check with Gemini now
            logging.warning(f'Check directly with GEMINI for {self.file_name}')
            em.get_obs_metadata(self._file_id)
            self._obs_id = ofr.repair_data_label(
                file_name, em.om.get('data_label'))
            clb = ofr.get_command_line_bits(self, self._obs_id, self.file_name)
            self._lineage = clb.lineage
            self._external_urls = clb.urls
        else:
            self._lineage = None
            self._external_urls = None
        if (self._file_id is None and self._obs_id is None and
                file_id is not None):
            self._file_id = file_id
            self._obs_id = file_id
        if file_id is not None:
            self._file_id = file_id
        logging.debug(self)

    def __str__(self):
        return f'obs_id {self._obs_id}, ' \
               f'file_id {self._file_id}, ' \
               f'file_name {self._file_name}, ' \
               f'lineage {self._lineage}, ' \
               f'external urls {self._external_urls}'

    def set_partial_args(self, pofr):
        temp = pofr.get_args(self._obs_id)
        if len(temp) == 1:
            self._lineage = temp[0].lineage
            self._external_urls = temp[0].urls
        else:
            raise mc.CadcException(
                f'Unexpected arguments for observation {self._obs_id} '
                f'file {self._file_id}')

    def _get_args(self):
        if self._lineage is None and self._external_urls is None:
            temp = em.get_gofr().get_args(self._obs_id)
            if len(temp) == 1:
                self._lineage = temp[0].lineage
                self._external_urls = temp[0].urls
                # format is GEMINI 'obs id'
                # use the repaired value
                self._obs_id = temp[0].obs_id.split()[1]
            else:
                if len(temp) == 0:
                    if self._file_id is None:
                        logging.warning(
                            f'Try a direct query for {self._obs_id}')
                        metadata = scrape.find_direct(self._obs_id)
                        if metadata is None or len(metadata) == 0:
                            raise mc.CadcException(
                                f'Direct query for obs id {self._obs_id}'
                                f' failed at Gemini')
                        else:
                            # assume filename from 0th entry, which
                            # should be ok, since the external urls and
                            # lineage will use all file names present
                            # in the returned metadata
                            self._file_id = GemName.remove_extensions(
                                metadata[0].get('filename'))
                            self._file_name = f'{self._file_id}.fits'
                            em.om.add(metadata, self._file_id)
                            temp_l = []
                            temp_e = []
                            for entry in metadata:
                                f_id = GemName.remove_extensions(
                                    entry.get('filename'))
                                f_name = f'{f_id}.fits'
                                temp_l.append(mc.get_lineage(
                                    ARCHIVE, f_id, f_name, SCHEME))
                                temp_e.append(f'{ofr.HEADER_URL}{f_name}')
                            self._lineage = ' '.join(ii for ii in temp_l)
                            self._external_urls = ' '.join(ii for ii in temp_e)
                            self._obs_id = ofr.repair_data_label(
                                self._file_name, self._obs_id)
                    else:
                        # Gemini obs id values are repaired from what
                        # archive.gemini.edu publishes, so check for the
                        # un-repaired value
                        repaired = em.get_gofr().get_obs_id(self._file_id)
                        # repaired = em.get_repaired_obs_id(self._file_id)
                        x = em.get_gofr().get_args(repaired)
                        self._lineage = x[0].lineage
                        self._external_urls = x[0].urls
                else:
                    found = False
                    for bits in temp:
                        urls = bits.urls.split()
                        for url in urls:
                            if self._file_name is None:
                                if (self._obs_id ==
                                        bits.obs_id.split()[1].strip()):
                                    logging.debug(f'Using existing obs id '
                                                  f'with {self._obs_id}')
                                    self._external_urls = bits.urls
                                    self._lineage = bits.lineage
                                    found = True
                                    break
                            elif url.endswith(self._file_name):
                                self._obs_id = bits.obs_id.split()[1]
                                logging.debug(f'Replaced obs id with '
                                              f'{self._obs_id}')
                                self._external_urls = bits.urls
                                self._lineage = bits.lineage
                                found = True
                                break
                    if not found:
                        raise mc.CadcException(
                            'Could not find obs id for file name {}'.format(
                                self._file_name))
        logging.debug(self)

    @property
    def file_uri(self):
        return '{}:{}/{}'.format(SCHEME, self.collection, self._file_name)

    @property
    def file_name(self):
        return self._file_name

    @file_name.setter
    def file_name(self, value):
        self._file_name = value.replace('.bz2', '').replace('.header', '')

    @property
    def compressed_file_name(self):
        return None

    @property
    def prev(self):
        return '{}.jpg'.format(self._file_id)

    @property
    def thumb(self):
        return '{}_th.jpg'.format(self._file_id)

    @property
    def file_id(self):
        return self._file_id

    @file_id.setter
    def file_id(self, value):
        self._file_id = value

    @property
    def lineage(self):
        if self._lineage is None:
            self._get_args()
        return self._lineage

    @property
    def external_urls(self):
        if self._external_urls is None:
            self._get_args()
        return self._external_urls

    @property
    def thumb_uri(self):
        """Note the 'ad' scheme - the thumbnail is generated at CADC,
        so acknowledge that with the ad URI."""
        return 'ad:{}/{}'.format(self.archive, self.thumb)

    def is_valid(self):
        return True

    @staticmethod
    def get_file_id(file_name):
        return GemName.remove_extensions(file_name)

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        # Note the .gz extension is on some TRECS files, not that it is
        # an accepted GEMINI extension
        return name.replace('.fits', '').replace('.bz2', ''). \
            replace('.header', '').replace('.jpg', '').replace('.gz', '')

    @staticmethod
    def is_preview(entry):
        return '.jpg' in entry


class GemNameBuilder(GemName):
    """Naming rules:
    - support mixed-case file name storage, exception for extensions, and
            mixed-case obs id values - the case the inputs are provided in are
            assumed to be correct.
    - support uncompressed files in storage
    """

    GEM_NAME_PATTERN = '*'

    def __init__(self,  obs_id, file_name, last_modified_s):
        # in this case there is no need to provide an init signature that's
        # consistent with other StorageName extensions, because the
        # constructor is used only in the builder class,
        # which is aware of Gemini data label/file name relationships
        logging.debug(f'parameters file_name {file_name} obs id {obs_id}')
        super(GemNameBuilder, self).__init__(
            obs_id=obs_id,fname_on_disk=file_name)
        # purposefully do not set self._file_name, since that indicates
        # later in processing that the file already exists in ad, and
        # naming patterns can be found by querying there
        self._file_id = GemName.remove_extensions(file_name)
        self._obs_id = ofr.repair_data_label(file_name, obs_id)
        self._last_modified_s = last_modified_s
        self._lineage = mc.get_lineage(
            ARCHIVE, self._file_id, file_name, SCHEME)
        self._external_urls = f'{ofr.HEADER_URL}{file_name}'
        self._file_name = None
        logging.debug(self)

    @property
    def external_urls(self):
        return self._external_urls

    @property
    def last_modified_s(self):
        return self._last_modified_s

    @property
    def lineage(self):
        return self._lineage
