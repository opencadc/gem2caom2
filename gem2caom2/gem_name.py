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
from gem2caom2.util import Inst, COLLECTION, SCHEME, V_SCHEME


__all__ = ['GemName']


HEADER_URL = 'https://archive.gemini.edu/fullheader/'


class GemName(mc.StorageName):
    """Naming rules:
    - support mixed-case file name storage, except for extensions, and
            mixed-case obs id values - the case the inputs are provided in are
            assumed to be correct.
    - support uncompressed files in storage
    - from Oliver Oberdorf at Gemini, 26-05-21: "prefer lower case, as that's
      what we use today"
      e.g.The observation ID GN-CAL20140811-1-033-BIAS has two artifacts:
            - N20140811S0033_BIAS.fits
            - N20140811S0033_bias.fits


    ALOPEKE/ZORRO::

    DB 31-08-20
    DATALAB can NOT be used for the CAOM2 Observation ID since it appears that
    the DATALAB value is identical for all files obtained for a single
    program. e.g. if the program ID is GN-2020A-DD-115 then the DATALAB value
    is always GN-2020A-DD-115-0-0.

    Instead, use the root of the filename as the observation ID.  e.g.
    N20200819A0003r.fits and N20200819A0003b.fits are two files generated from
    a single observation (r = red channel, b = blue channel).  Use
    N20200819A0003 as the observation ID with two planes given by the two
    colours of data.

    DB 01-09-20
    Gemini has kludged the headers so that every observation for a single
    program has the same DATALAB in the header.  This is what we usually use
    for the observation ID.  Each single ‘observation’ actually produces two
    files (not a single MEF file) for the red and blue channels so to me it
    would make the most sense to group these two files as a single observation
    with two artifacts given by uri’s pointing to the two files.  And this is
    a single plane, correct?

    PD 01-09-20
    What is the meaning of red and blue channels? different energy bands?

    DB 02-09-20
    Yes.  there’s a dichroic that directs light shortward of 675nm to one
    detector (through one of several possible filters) and light longward of
    675nm to a second detector (through another filter).   But instead of
    generating a single MEF file they generate two files, e.g.
    N20191219A0004b.fits and N20191219A0004r.fits.

    PD 02-09-20
    This seems very much like MACHO... if those two files are images in the
    normal sense then it could make sense to create separate planes with
    dataProductType = image that end up with the correct (distinct) energy
    metadata. It is OK for an observation to create two sibling products and
    two planes probably captures the goal of this instrument/observing mode
    more directly.

    PD 15-12-20 IPM
    When working with the new storage system at CADC:
    Artifact URIs for files, previews, obtained from archive.gemini.edu
    should be:
      gemini:GEMINI/file_name."fits|jpg"[.bz2]
    Artifact URIs for thumbnails created at CADC should be:
      cadc:GEMINI/file_name_th.jpg
    """

    GEM_NAME_PATTERN = '*'

    def __init__(
        self,
        file_name=None,
        obs_id=None,
        instrument=None,
        entry=None,
    ):
        if instrument in [Inst.ALOPEKE, Inst.ZORRO]:
            self._file_name = file_name
            self._file_id = GemName.remove_extensions(self._file_name)
            self._obs_id = self._file_id[:-1]
            self._product_id = self._file_id
            super(GemName, self).__init__(
                obs_id=self._obs_id,
                collection=COLLECTION,
                collection_pattern=GemName.GEM_NAME_PATTERN,
                fname_on_disk=self.file_name,
                scheme=SCHEME,
                entry=entry,
            )
        else:
            # try to set the file name, if that information is available

            # file_name is assumed to be the file name in ad
            # because the GEM files are stored uncompressed,
            # while the files available from Gemini are bz2.
            self._file_name = None
            self._file_id = None
            if file_name is not None:
                self._file_id = GemName.get_file_id(file_name)
                self.file_name = file_name
            if obs_id is not None:
                self._obs_id = obs_id
            super(GemName, self).__init__(
                obs_id=obs_id,
                collection=COLLECTION,
                collection_pattern=GemName.GEM_NAME_PATTERN,
                fname_on_disk=self.file_name,
                scheme=SCHEME,
                entry=entry,
            )
            if self._obs_id is None:
                temp = em.get_gofr().get_obs_id(self._file_id)
                if temp is not None:
                    self._obs_id = GemName.remove_extensions(temp)
            if (
                self._fname_on_disk is None
                and self._file_name is None
                and self._obs_id is None
            ):
                raise mc.CadcException('Require a name.')
            self._product_id = self._file_id
        self._v_scheme = V_SCHEME
        self._v_collection = COLLECTION
        self._logger = logging.getLogger(__name__)
        self._logger.debug(self)

    def __str__(self):
        return (
            f'\n'
            f'     obs_id:{self._obs_id}\n'
            f'    file_id:{self._file_id}\n'
            f'   file_name:{self._file_name}\n'
            f'    file_uri:{self.file_uri}\n'
            f'   thumb_uri:{self.thumb_uri}\n'
            f'    prev_uri:{self.prev_uri}\n'
            f'source_names:{self._source_names}\n'
        )

    @property
    def file_uri(self):
        return f'{self.scheme}:{self._v_collection}/{self._file_name}'

    @property
    def file_name(self):
        return self._file_name

    @file_name.setter
    def file_name(self, value):
        if value is None:
            self._file_name = None
        else:
            self._file_name = value.replace('.bz2', '').replace('.header', '')

    @property
    def compressed_file_name(self):
        return None

    @property
    def prev(self):
        return f'{self._file_id}.jpg'

    @property
    def thumb(self):
        return f'{self._file_id}_th.jpg'

    @property
    def file_id(self):
        return self._file_id

    @file_id.setter
    def file_id(self, value):
        self._file_id = value

    @property
    def lineage(self):
        if '_th.jpg' in self._file_name:
            # thumbnail
            return mc.get_lineage(
                self._v_collection,
                self.product_id,
                self._file_name,
                self._v_scheme,
            )
        else:
            return mc.get_lineage(
                self._v_collection,
                self.product_id,
                self._file_name,
                SCHEME,
            )

    @property
    def external_urls(self):
        return f'{HEADER_URL}{self._file_id}.fits'

    @property
    def prev_uri(self):
        return f'{self.scheme}:{self._v_collection}/{self.prev}'

    @property
    def product_id(self):
        return self._product_id

    @property
    def product_id(self):
        return self._file_id

    @property
    def thumb_uri(self):
        """Note the v_scheme - the thumbnail is generated at CADC,
        so acknowledge that with the CADC URI."""
        return f'{self._v_scheme}:{self._v_collection}/{self.thumb}'

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
        return (
            name.replace('.fits', '')
            .replace('.bz2', '')
            .replace('.header', '')
            .replace('.jpg', '')
            .replace('.gz', '')
        )

    @staticmethod
    def is_preview(entry):
        return '.jpg' in entry
