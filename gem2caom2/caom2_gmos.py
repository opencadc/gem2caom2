# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
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

import re
import caom2
from gem2caom2.svofps import filter_metadata


# __all__ = ['GMOS']


class GMOS():
    """
    GMOS-N/S CAOM2 metadata generation
    """

    def __init__(self):
        """
        GMOS N/S resolving power, dispersion, etc. from Gemini web pages
        """
        self.energy_band = caom2.EnergyBand['OPTICAL']

        self.res_power = {
            'B1200': 3744.0,
            'R831': 4396.0,
            'B600': 1688.0,
            'R600': 3744.0,
            'R400': 1918.0,
            'R150': 631.0
        }

        # Angstroms/pixel
        self.dispersion = {
            'B1200': 0.245,
            'R831': 0.36,
            'B600': 0.475,
            'R600': 0.495,
            'R400': 0.705,
            'R150': 1.835
        }

    def energy_metadata(self, obs_metadata):
        energy_metadata = {
            'energy': True,
            'energy_band': self.energy_band
        }

        # Determine energy metadata for the plane.  
        # No energy information is determined for biases or darks.  The
        # latter are sometimes only identified by a 'blank' filter.  e.g.
        # NIRI 'flats' are sometimes obtained with the filter wheel blocked off.

        if obs_metadata['observation_type'] in ('BIAS', 'DARK'):
            energy_metadata['energy'] = False
            return energy_metadata

        if obs_metadata['mode'] == 'imaging':
            filter_md = filter_metadata(obs_metadata['instrument'],
                                        obs_metadata['filter_name'])

            delta = filter_md['wl_eff_width']
            reference_wavelength = filter_md['wl_eff']
            resolving_power = reference_wavelength/delta
            reference_wavelength /= 1.0e10
            delta /= 1.0e10
        elif obs_metadata['mode'] in ('LS', 'spectroscopy', 'MOS'):
            reference_wavelength = obs_metadata['central_wavelength']
            # Ignore energy information if value of 'central_wavelength' = 0.0
            if reference_wavelength == 0.0:
                energy_metadata['energy'] = False
                # if 'focus' in fpmask:
                #     md['observation_type'] = 'FOCUS'
                return energy_metadata

            resolving_power = self.res_power[obs_metadata['disperser']]
            delta = self.dispersion[obs_metadata['disperser']]

            reference_wavelength /= 1.0e6
            delta /= 1.0e10

        filter_name = re.sub(r'&', ' & ', obs_metadata['filter_name'])
        energy_metadata['filter_name'] = filter_name
        energy_metadata['wavelength_type'] = 'WAVE'
        energy_metadata['wavelength_unit'] = 'm'
        energy_metadata['number_pixels'] = 1
        energy_metadata['reference_wavelength'] = reference_wavelength
        energy_metadata['delta'] = delta
        energy_metadata['resolving_power'] = resolving_power
        energy_metadata['reference_pixel'] = 1.0

        # On occasion (e.g. Moon observations) two filters with
        # inconsistent passbands result in negative r.  e.g. RG610 + g
        # Skip energy metadata in this case.
        if energy_metadata['resolving_power'] < 0:
            energy_metadata['energy'] = False
  
        return energy_metadata
