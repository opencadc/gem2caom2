# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2024.                            (c) 2024.
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
#  Revision: 4
#
# ***********************************************************************
#

from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from caom2 import ReleaseType
from caom2pipe.manage_composable import PreviewVisitor
from gem2caom2.util import Inst
from matplotlib import colors as colors
from os.path import basename

import matplotlib.pyplot as plt
import numpy as np


class GHOSTPreviews(PreviewVisitor):

    def __init__(self, **kwargs):
        super().__init__(ReleaseType.META, **kwargs)

    def generate_plots(self, obs_id):
        self._logger.debug(f'Begin generate_plots for {obs_id}')
        # code by https://github.com/dbohlender
        plt.style.use(astropy_mpl_style)
        fig = plt.figure(figsize=(8, 8))
        hdulist = fits.open(self._science_fqn)
        target = hdulist[0].header['OBJECT']
        # The extension numbers must be determined from the extension with FITS header values for NAXIS = 2 and
        # CAMERA = RED|BLUE. Each channel normally has 4 image exiensions.  If multiple exposures are stored in
        # the file then only the first is used for the preview.
        blue_ext = []
        red_ext = []
        # Run through all of the extensions to fine the RED and BLUE extensions.
        for h in range(len(hdulist)):
            naxis = hdulist[h].header['NAXIS']
            if naxis == 2:
                camera = hdulist[h].header['CAMERA']
                if camera == 'RED':
                    red_ext.append(h)
                if camera == 'BLUE':
                    blue_ext.append(h)
        red_data = []
        blue_data = []
        for i in red_ext:
            red_data.append(fits.getdata(self._science_fqn, ext=i))
        for i in blue_ext:
            blue_data.append(fits.getdata(self._science_fqn, ext=i))
        # Not very elegant, but the image arrays need to be combined to create a final large, 2D image for both
        # channels. They are then 'flip'ed appropriately so that short wavelengths are at the top of each image
        # Note:  this assumes that the order of each image extension does not change!
        if red_data and blue_data:
            red_image1 = np.concatenate((red_data[0], red_data[1]), axis=1)
            red_image2 = np.concatenate((red_data[3], red_data[2]), axis=1)
            red_image = np.concatenate((red_image1, red_image2), axis=0)
            red_image = np.flip(red_image, axis=1)
            blue_image1 = np.concatenate((blue_data[0], blue_data[1]), axis=1)
            blue_image2 = np.concatenate((blue_data[3], blue_data[2]), axis=1)
            blue_image = np.concatenate((blue_image1, blue_image2), axis=0)
            blue_image = np.flip(blue_image)
            fig.add_subplot(2, 1, 1)
            plt.axis('off')
            plt.title(f'{basename(self._science_fqn)}:   {target}\nBlue Channel')
            plt.imshow(blue_image, cmap='Blues_r', norm=colors.LogNorm())
            fig.add_subplot(2, 1, 2)
            plt.axis('off')
            plt.title(f'Red Channel')
            plt.imshow(red_image, cmap='Reds_r', norm=colors.LogNorm())
            plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=0.92, wspace=0.0, hspace=0.1)
            plt.savefig(self._preview_fqn)
            plt.close()
            self._logger.debug('Finish generate_plots')
            return self._save_figure()
        else:
            self._logger.warning(f'Found no image metadata for {self._storage_name.file_uri}')
            return 0


def visit(observation, **kwargs):
    if observation.instrument.name == Inst.GHOST.value:
        return GHOSTPreviews(**kwargs).visit(observation)
    else:
        return observation
