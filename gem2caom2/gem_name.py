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


__all__ = ['GemName', 'COLLECTION', 'ARCHIVE', 'SCHEME']


COLLECTION = 'GEMINI'
ARCHIVE = 'GEM'
SCHEME = 'gemini'
HEADER_URL = 'https://archive.gemini.edu/fullheader/'


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
        super(GemName, self).__init__(
            obs_id=obs_id, collection=ARCHIVE,
            collection_pattern=GemName.GEM_NAME_PATTERN,
            fname_on_disk=self.file_name,
            scheme=SCHEME)
        if self._obs_id is None:
            temp = em.get_gofr().get_obs_id(self._file_id)
            if temp is not None:
                self._obs_id = GemName.remove_extensions(temp)
        if (self._fname_on_disk is None and self._file_name is None and
                self._obs_id is None):
            raise mc.CadcException('Require a name.')
        if (self._file_id is None and self._obs_id is None and
                file_id is not None):
            self._file_id = file_id
            self._obs_id = file_id
        if file_id is not None:
            self._file_id = file_id
        self._logger = logging.getLogger(__name__)
        self._logger.debug(self)

    def __str__(self):
        return f'obs_id {self._obs_id}, ' \
               f'file_id {self._file_id}, ' \
               f'file_name {self._file_name}'

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
        return mc.get_lineage(ARCHIVE, self._file_id, self._file_name, SCHEME)

    @property
    def external_urls(self):
        return f'{HEADER_URL}{self._file_id}.fits'

    @property
    def product_id(self):
        return self._file_id

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
    def get_file_name_from(file_id):
        # makes the assumption that stored data is un-compressed
        return f'{file_id}.fits'

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
