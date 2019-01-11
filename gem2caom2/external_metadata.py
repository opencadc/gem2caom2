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

import os
import re
import requests
from requests.adapters import HTTPAdapter
from urllib3 import Retry

import caom2
from caom2pipe import manage_composable as mc
from gem2caom2.svofps import filter_metadata

GEMINI_METADATA_URL = 'https://archive.gemini.edu/jsonsummary/canonical/'
GEMINI_FITS_HEADER_URL = 'https://archive.gemini.edu/fullheader/'

GMOS_ENERGY_BAND = caom2.EnergyBand['OPTICAL']
NIRI_ENERGY_BAND = caom2.EnergyBand['INFRARED']

GMOS_RESOLVING_POWER = {
    'B1200': 3744.0,
    'R831': 4396.0,
    'B600': 1688.0,
    'R600': 3744.0,
    'R400': 1918.0,
    'R150': 631.0
}

# Angstroms/pixel
GMOS_DISPERSION = {
    'B1200': 0.245,
    'R831': 0.36,
    'B600': 0.475,
    'R600': 0.495,
    'R400': 0.705,
    'R150': 1.835
}

NIRI_RESOLVING_POWER = {
    'J': {},
    'H': {},
    'L': {},
    'M': {},
    'K': {
        'f6-2pix': 1300.0,
        'f6-4pix': 780.0,
        'f6-6pix': 520.0
    }
}


def get_fits_headers(file_name):
    """
    Get the headers for the given FITS file name.

    :param file_name: The FITS file name.
    :return: List of FITS headers.
    """
    # file_name should end in .fits, strip off extensions after .fits
    if not file_name.endswith('.fits'):
        file_name = os.path.splitext(file_name)[0]

    gemini_url = '{}{}'.format(GEMINI_FITS_HEADER_URL, file_name)

    # Open the URL and fetch the FITS headers for the observation
    session = requests.Session()
    retries = 10
    retry = Retry(total=retries, read=retries, connect=retries,
                  backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    try:
        response = session.get(gemini_url, timeout=20)
        header = response.text.split('\n')
        response.close()
    except Exception as e:
        raise mc.CadcException(
            'Unable to download Gemini observation header from {} because {}'
                .format(gemini_url, str(e)))
    return header


def gmos_metadata(obs_metadata):
    """
    Calculate GMOS energy metadata using the Gemini observation metadata.

    :param obs_metadata: Dictionary of observation metadata.
    :return: Dictionary of energy metadata
    """
    metadata = {
        'energy': True,
        'energy_band': GMOS_ENERGY_BAND
    }

    # Determine energy metadata for the plane.
    # No energy information is determined for biases or darks.  The
    # latter are sometimes only identified by a 'blank' filter.  e.g.
    # NIRI 'flats' are sometimes obtained with the filter wheel blocked off.
    if obs_metadata['observation_type'] in ('BIAS', 'DARK'):
        metadata['energy'] = False
        return metadata

    reference_wavelength = 0.0
    delta = 0.0
    resolving_power = 0.0

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
            metadata['energy'] = False
            return metadata

        resolving_power = GMOS_RESOLVING_POWER[obs_metadata['disperser']]
        delta = GMOS_DISPERSION[obs_metadata['disperser']]

        reference_wavelength /= 1.0e6
        delta /= 1.0e10

    filter_name = re.sub(r'&', ' & ', obs_metadata['filter_name'])
    metadata['filter_name'] = filter_name
    metadata['wavelength_type'] = 'WAVE'
    # metadata['wavelength_unit'] = 'm'
    metadata['number_pixels'] = 1
    metadata['reference_wavelength'] = reference_wavelength
    metadata['delta'] = delta
    metadata['resolving_power'] = resolving_power
    metadata['reference_pixel'] = 1.0

    # On occasion (e.g. Moon observations) two filters with
    # inconsistent passbands result in negative r.  e.g. RG610 + g
    # Skip energy metadata in this case.
    if metadata['resolving_power'] < 0:
        metadata['energy'] = False

    return metadata


def niri_metadata(obs_metadata):
    """
    Calculate NIRI energy metadata using the Gemini observation metadata.

    :param obs_metadata: Dictionary of observation metadata.
    :return: Dictionary of energy metadata
    """
    metadata = {
        'energy': True,
        'energy_band': NIRI_ENERGY_BAND
    }

    # Determine energy metadata for the plane.
    # No energy information is determined for darks.  The
    # latter are sometimes only identified by a 'blank' filter.  e.g.
    # NIRI 'flats' are sometimes obtained with the filter wheel blocked off.
    headers = get_fits_headers(obs_metadata['filename'])
    header_filters = []
    filters2ignore = ['open', 'INVALID', 'PK50', 'pupil']
    for header in headers:
        if 'FILTER' in header:
            if any(x in header for x in filters2ignore):
                continue
            else:
                filter = "".join(re.findall(r'\'(.+?)\'', header))
                filter = filter.replace('_', '-').strip()
                filter = ''.join('' if ch in '()' else ch for ch in filter)
                header_filters.append(filter)
        filters = "&".join(header_filters)

    if obs_metadata['observation_type'] in 'DARK' or 'blank' in filters:
        metadata['energy'] = False
        return metadata

    reference_wavelength = 0.0
    delta = 0.0
    resolving_power = 0.0

    filter_md = filter_metadata(obs_metadata['instrument'], filters)
    if obs_metadata['mode'] == 'imaging':

        delta = filter_md['wl_eff_width']
        reference_wavelength = filter_md['wl_eff']
        resolving_power = reference_wavelength/delta
        reference_wavelength /= 1.0e10
        delta /= 1.0e10
    # elif obs_metadata['mode'] in ('LS', 'spectroscopy'):
    #    # this code has to be rewritten for NIRI!!!
    #    reference_wavelength = obs_metadata['central_wavelength']
    #    nrgdim = int(niri_metadata['naxis2']/bin_y)

    #    # Ignore energy information if value of 'central_wavelength' = 0.0
    #    if reference_wavelength == 0.0 \
    #            or obs_metadata['observation_type'] == 'BIAS':
    #        metadata['energy'] = False
    #        return metadata

    #   if 'focus' in fpmask:
    #       obstype = 'FOCUS'
    #    resolving_power = NIRI_RESOLVING_POWER[bandpassname][fpmask]
    #    delta = filter_md['wl_eff_width']/metadata['naxis1']
    #    reference_wavelength /= 1.0e6
    #    delta /= 1.0e10

    filters = re.sub(r'&', ' & ', filters)
    filters = re.sub(r'-G.{4}(|w)', '', filters)
    metadata['filter_name'] = filters
    metadata['wavelength_type'] = 'WAVE'
    # metadata['wavelength_unit'] = 'm'
    metadata['number_pixels'] = 1024
    metadata['reference_wavelength'] = reference_wavelength
    metadata['delta'] = delta
    metadata['resolving_power'] = resolving_power
    metadata['reference_pixel'] = 512.0

    return metadata
