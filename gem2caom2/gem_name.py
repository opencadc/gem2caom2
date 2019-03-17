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
from caom2pipe import manage_composable as mc
from caom2pipe import execute_composable as ec


__all__ = ['GemName', 'COLLECTION', 'ARCHIVE', 'SCHEME']


COLLECTION = 'GEMINI'
ARCHIVE = 'GEM'
SCHEME = 'gemini'


class GemName(ec.StorageName):
    """Naming rules:
    - support mixed-case file name storage, exception for extensions, and
            mixed-case obs id values - the case the inputs are provided in are
            assumed to be correct.
    - support uncompressed files in storage
    """

    GEM_NAME_PATTERN = '*'

    def __init__(self, fname_on_disk=None, file_name=None, obs_id=None):
        if file_name is not None:
            self.file_id = GemName.get_file_id(file_name)
            if '.fits' in file_name:
                self.fname_in_ad = '{}.fits'.format(self.file_id)
            elif GemName.is_preview(file_name):
                self.fname_in_ad = '{}.jpg'.format(self.file_id)
            else:
                raise mc.CadcException(
                    'Unrecognized file name format {}'.format(file_name))
        elif fname_on_disk is not None:
            self.file_id = GemName.get_file_id(fname_on_disk)
            if '.fits' in fname_on_disk:
                self.fname_in_ad = '{}.fits'.format(self.file_id)
            elif GemName.is_preview(fname_on_disk):
                self.fname_in_ad = '{}.jpg'.format(self.file_id)
            else:
                raise mc.CadcException(
                    'Unrecognized file name format {}'.format(fname_on_disk))
        else:
            raise mc.CadcException('Require file name.')
        super(GemName, self).__init__(
            obs_id=None, collection=ARCHIVE,
            collection_pattern=GemName.GEM_NAME_PATTERN,
            fname_on_disk=fname_on_disk,
            scheme=SCHEME)
        self.obs_id = obs_id

    @property
    def file_uri(self):
        return '{}:{}/{}'.format(SCHEME, self.collection, self.file_name)

    @property
    def file_name(self):
        return self.fname_in_ad

    @property
    def compressed_file_name(self):
        return None

    @property
    def prev(self):
        return '{}.jpg'.format(GemName.get_file_id(self.fname_in_ad))

    @property
    def thumb(self):
        return '{}_th.jpg'.format(GemName.get_file_id(self.fname_in_ad))

    @property
    def obs_id(self):
        return self._obs_id

    @obs_id.setter
    def obs_id(self, value):
        self._obs_id = value

    @property
    def file_id(self):
        return self._file_id

    @file_id.setter
    def file_id(self, value):
        self._file_id = value

    def is_valid(self):
        return True

    # def _get_obs_id(self):
    #     if 'TMP' in self.file_uri:
    #         # TODO for testing only
    #         with open(self.file_uri) as f:
    #             headers = f.readlines()
    #     elif '.fits' in self.file_uri:
    #         headers = mc.get_cadc_headers(self.file_uri)
    #     else:
    #         temp_uri = self.file_uri.replace('.jpg', '.fits')
    #         headers = mc.get_cadc_headers(temp_uri)
    #     fits_headers = ac.make_headers_from_string(headers)
    #     obs_id = fits_headers[0].get('DATALAB')
    #     return obs_id

    @staticmethod
    def get_file_id(file_name):
        return GemName.remove_extensions(file_name)

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        return name.replace('.fits', '').replace('.bz2', ''). \
            replace('.header', '').replace('.jpg', '')

    @staticmethod
    def is_preview(entry):
        return '.jpg' in entry
