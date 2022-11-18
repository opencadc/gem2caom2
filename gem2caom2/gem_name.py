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

from gem2caom2.obs_file_relationship import remove_extensions


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
    ):
        super().__init__(file_name=file_name.replace('.header', ''))
        # use the file id because the extension doesn't help much in the archive.gemini.edu URL
        self._source_names = [self._file_id]

    @property
    def prev(self):
        return f'{self._file_id}.jpg'

    @property
    def thumb(self):
        return f'{self._file_id}_th.jpg'

    @property
    def prev_uri(self):
        # use the 'gemini' scheme because the previews are from archive.gemini.edu
        return self._get_uri(self.prev, mc.StorageName.scheme)

    @property
    def thumb_uri(self):
        # use the 'cadc' scheme because the thumbnails are generated at CADC from the archive.gemini.edu previews
        return self._get_uri(self.thumb, mc.StorageName.preview_scheme)

    def is_valid(self):
        return True

    def set_file_id(self):
        self._file_id = remove_extensions(self._file_name)

    def set_obs_id(self):
        if self._file_id[-1] in ['b', 'r']:
            self._obs_id = self._file_id[:-1]

    def set_product_id(self):
        if self._file_id[-1] in ['b', 'r']:
            self._product_id = self._file_id
        else:
            if self._file_id.startswith('SDC'):
                # DB 20-07-21
                #  each pair of H/K files will be one observation with one
                #  plane with two artifacts.
                self._product_id = self._file_id.replace(
                    'SDCH', 'SDC'
                ).replace('SDCK', 'SDC')
            else:
                self._product_id = self._file_id

    def set_destination_uris(self):
        self._destination_uris = [self.file_uri]

    @staticmethod
    def get_file_name_from(file_id):
        # makes the assumption that stored data is un-compressed
        return f'{file_id}.fits'

    @staticmethod
    def is_preview(entry):
        return '.jpg' in entry
