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

from gem2caom2 import external_metadata as em
from gem2caom2.obs_file_relationship import GemObsFileRelationship


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
        if em.gofr is None:
            em.gofr = GemObsFileRelationship('/app/data/from_paul.txt')

        # try to set the file name, if that information is available

        # file_name is assumed to be the file name in ad
        # because the GEM files are stored uncompressed,
        # while the files available from Gemini are bz2.
        self._file_name = None
        if file_name is not None:
            self._file_id = GemName.get_file_id(file_name)
            self.file_name = file_name
        elif fname_on_disk is not None:
            self._file_id = GemName.get_file_id(fname_on_disk)
            self.file_name = fname_on_disk
        elif obs_id is not None:
            self._obs_id = obs_id
        else:
            raise mc.CadcException('Require file name.')
        super(GemName, self).__init__(
            obs_id=obs_id, collection=ARCHIVE,
            collection_pattern=GemName.GEM_NAME_PATTERN,
            fname_on_disk=fname_on_disk,
            scheme=SCHEME)
        self.gofr_file_names = None
        self.gofr_clb = None

    @property
    def file_uri(self):
        return '{}:{}/{}'.format(SCHEME, self.collection, self._file_name)

    @property
    def file_name(self):
        return self._file_name

    @file_name.setter
    def file_name(self,value):
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
        if self.gofr_clb is None:
            self.gofr_clb = em.gofr.get_args(self._obs_id)
        if len(self.gofr_clb) > 1:
            raise mc.CadcException(
                'Too many obs ids for {}'.format(self._obs_id))
        return self.gofr_clb[0].lineage

    @property
    def external_urls(self):
        if self.gofr_clb is None:
            self.gofr_clb = em.gofr.get_args(self._obs_id)
        if len(self.gofr_clb) > 1:
            raise mc.CadcException(
                'Too many obs ids for {}'.format(self._obs_id))
        return self.gofr_clb[0].urls

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
        return name.replace('.fits', '').replace('.bz2', ''). \
            replace('.header', '').replace('.jpg', '')

    @staticmethod
    def is_preview(entry):
        return '.jpg' in entry
