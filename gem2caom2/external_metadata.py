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

# import os
import logging
import re
import requests
from requests.adapters import HTTPAdapter
from urllib3 import Retry

from bs4 import BeautifulSoup

import caom2
from caom2pipe import manage_composable as mc
from gem2caom2.svofps import filter_metadata
from gem2caom2 import gemini_obs_metadata as gom


GEMINI_METADATA_URL = \
    'https://archive.gemini.edu/jsonsummary/canonical/filepre='
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

# values from
# https://www.gemini.edu/sciops/instruments/niri/spectroscopy/grisms
NIRI_RESOLVING_POWER = {
    'J': {
        'f6-2pix': 770.0,
        'f6-4pix': 610.0,
        'f6-6pix': 460.0,
        'f6-2pixB1': 770.0,
        'f6-4pixB1': 650.0,
        'f6-6pixB1': 480.0,
        'f32-4pix': 1000.0,
        'f32-6pix': 620.0,  # f32-7pix
        'f32-9pix': 450.0  # f32-10pix
    },
    'H': {
        'f6-2pix': 1650.0,
        'f6-4pix': 825.0,
        'f6-6pix': 520.0,
        'f6-2pixB1': 1650.0,
        'f6-4pixB1': 940.0,
        'f6-6pixB1': 550.0,
        'f32-4pix': 880.0,
        'f32-6pix': 630.0,  # f32-7pix
        'f32-9pix': 500.0  # f32-10pix
    },
    'L': {
        'f6-2pix': 1100.0,
        'f6-4pix': 690.0,
        'f6-6pix': 460.0,
        'f6-2pixB1': 1100.0,
        'f6-4pixB1': 770.0,
        'f6-6pixB1': 490.0,
    },
    'M': {
        'f6-2pix': 1100.0,
        'f6-4pix': 770.0,
        'f6-6pix': 460.0
    },
    'K': {
        'f6-2pix': 1300.0,
        'f6-4pix': 780.0,
        'f6-6pix': 520.0,
        'f32-4pix': 1280.0,
        'f32-6pix': 775.0,  # f32-7pix
        'f32-9pix': 570.0  # f32-10pix
    }
}

# select filter_id, wavelength_central, wavelength_lower, wavelength_upper
# from gsa..gsa_filters where instrument = 'NICI'
# 0 - central
# 1 - lower
# 2 - upper

PHOENIX = {'2030 (4)': [4.929000, 4.808000, 5.050000],
           '2030 (9)': [4.929000, 4.808000, 5.050000],
           '2150 (2)': [4.658500, 4.566000, 4.751000],
           '2462 (5)': [4.078500, 4.008000, 4.149000],
           '2734 (4)': [3.670500, 3.610000, 3.731000],
           '2870 (7)': [3.490500, 3.436000, 3.545000],
           '3010 (4)': [3.334500, 3.279000, 3.390000],
           '3100 (11)': [3.240000, 3.180000, 3.300000],
           '3290 (7)': [3.032500, 2.980000, 3.085000],
           '4220 (5)': [2.370000, 2.348000, 2.392000],
           '4308 (6)': [2.322500, 2.296000, 2.349000],
           '4396 (8)': [2.272500, 2.249000, 2.296000],
           '4484 (8)': [2.230000, 2.210000, 2.250000],
           '4578 (6)': [2.185000, 2.160000, 2.210000],
           '4667 (9)': [2.143000, 2.120000, 2.166000],
           '4748 (11)': [2.104000, 2.082000, 2.126000],
           '6073 (10)': [1.647000, 1.632000, 1.662000],
           '6420 (12)': [1.557500, 1.547000, 1.568000],
           '7799 (10)': [1.280500, 1.271000, 1.290000],
           '8265 (13)': [1.204500, 1.196000, 1.213000],
           '9232 (3)': [1.083000, 1.077000, 1.089000],
           'L2870 (7)': [3.490500, 3.436000, 3.545000]}

obs_metadata = {}
om = None
fm = {}


def get_obs_metadata(file_id):
    """
    Download the Gemini observation metadata for the given obs_id.

    :param file_id: The file ID
    :return: Dictionary of observation metadata.
    """
    logging.debug('Begin get_obs_metadata')
    gemini_url = '{}{}'.format(GEMINI_METADATA_URL, file_id)

    # Open the URL and fetch the JSON document for the observation
    session = requests.Session()
    retries = 10
    retry = Retry(total=retries, read=retries, connect=retries,
                  backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    try:
        response = session.get(gemini_url, timeout=20)
        metadata = response.json()[0]
        response.close()
    except Exception as e:
        raise mc.CadcException(
            'Unable to download Gemini observation metadata from {} because {}'
                .format(gemini_url, str(e)))
    global obs_metadata
    obs_metadata = metadata
    global om
    om = gom.GeminiObsMetadata(metadata, file_id)
    logging.debug('End get_obs_metadata')


def get_pi_metadata(program_id):
    program_url = 'https://archive.gemini.edu/programinfo/' + program_id

    # Open the URL and fetch the JSON document for the observation
    session = requests.Session()
    retries = 10
    retry = Retry(total=retries, read=retries, connect=retries,
                  backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    try:
        response = session.get(program_url, timeout=20)
        xml_metadata = response.text
        response.close()
        logging.error('got pi metdata')
    except Exception as e:
        raise mc.CadcException(
            'Unable to download Gemini observation metadata from {} because {}'
                .format(program_url, str(e)))
    metadata = None
    soup = BeautifulSoup(xml_metadata, 'lxml')
    tds = soup.find_all('td')
    if len(tds) > 0:
        title = tds[1].contents[0].replace('\n', ' ')
        pi_name = tds[3].contents[0]
        metadata = {'title': title,
                    'pi_name': pi_name}
    logging.debug('End get_obs_metadata')
    return metadata


def get_filter_metadata(instrument, filter_name):
    """A way to lazily initialize all the filter metadata reads from SVO."""
    global fm
    if instrument in fm and filter_name in fm[instrument]:
        result = fm[instrument][filter_name]
    else:
        result = filter_metadata(instrument, filter_name)
        if instrument in fm:
            temp = fm[instrument]
            temp[filter_name] = result
        else:
            fm[instrument] = {filter_name: result}
    return result


def gmos_metadata():
    """
    Calculate GMOS energy metadata using the Gemini observation metadata.
    Imaging observations require a filter lookup to an external service.

    :param obs_metadata: Dictionary of observation metadata.
    :return: Dictionary of energy metadata
    """
    logging.debug('Begin gmos_metadata')
    metadata = {
        'energy': True,
        'energy_band': GMOS_ENERGY_BAND
    }

    # Determine energy metadata for the plane.
    # No energy information is determined for biases or darks.  The
    # latter are sometimes only identified by a 'blank' filter.  e.g.
    # NIRI 'flats' are sometimes obtained with the filter wheel blocked off.
    global obs_metadata
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

    filter_name = re.sub(r'\+', ' + ', obs_metadata['filter_name'])
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

    logging.debug('End gmos_metadata')
    return metadata
