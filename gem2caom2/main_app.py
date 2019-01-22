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
"""
Notes on the GEM archive/GEMINI collection:

1. Must use the file name as the starting point for work, because that's
what is coming back from the ad query, and the ad query is what is being
used to trigger the work.

2. Must find the observation ID value from the file header information,
because the observation ID is how to get the existing CAOM instance.

3. Artifact URIs in existing observations reference the gemini schema, not
the ad schema.

4. TODO - what happens to obtaining the observation ID, if a preview
file is retrieved prior to the FITS file being retrieved?

Because of this, make the GemName a class that, standing on it's own,
can retrieve the observation ID value from the headers for a file.

"""
import importlib
import logging
import os
import sys
import traceback
import requests
from requests.adapters import HTTPAdapter
from urllib3 import Retry

from caom2 import Observation, ObservationIntentType, DataProductType
from caom2 import CalibrationLevel, TargetType, ProductType
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import manage_composable as mc
from caom2pipe import execute_composable as ec
from caom2pipe import astro_composable as ac

from gem2caom2.external_metadata import gmos_metadata, niri_metadata


__all__ = ['main_app2', 'update', 'GemName', 'COLLECTION', 'APPLICATION',
           'SCHEME', 'ARCHIVE']


APPLICATION = 'gem2caom2'
COLLECTION = 'GEMINI'
ARCHIVE = 'GEM'
SCHEME = 'gemini'
GEMINI_METADATA_URL = 'https://archive.gemini.edu/jsonsummary/canonical/'

obs_metadata = {}


class GemName(ec.StorageName):
    """Naming rules:
    - support upper-case file name storage, exception for extensions, and
            upper-case obs id values
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
        return GemName.remove_extensions(file_name.lower()).upper()

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        return name.replace('.fits', '').replace('.gz', ''). \
            replace('.header', '').replace('.jpg', '')

    @staticmethod
    def is_preview(entry):
        return '.jpg' in entry


def get_obs_metadata(obs_id):
    """
    Download the Gemini observation metadata for the given obs_id.

    :param obs_id: The Obs ID
    :return: Dictionary of observation metadata.
    """
    gemini_url = '{}{}'.format(GEMINI_METADATA_URL, obs_id)

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
        logging.error('got obs metdata')
    except Exception as e:
        raise mc.CadcException(
            'Unable to download Gemini observation metadata from {} because {}'
                .format(gemini_url, str(e)))
    return metadata


def get_energy_metadata():
    """
    For the given observation retrieve the energy metadata.

    :return: Dictionary of energy metadata.
    """
    logging.debug('Begin get_energy_metadata')
    global obs_metadata
    instrument = obs_metadata['instrument']
    if instrument in ['GMOS-N', 'GMOS-S']:
        energy_metadata = gmos_metadata(obs_metadata)
    elif instrument in ['NIRI']:
        energy_metadata = niri_metadata(obs_metadata)
    elif instrument in ['GNIRS']:
        # TODO
        energy_metadata = niri_metadata(obs_metadata)
    else:
        raise mc.CadcException(
            'Do not understand energy for instrument {}'.format(instrument))
    logging.debug(
        'End get_energy_metadata for instrument {}'.format(instrument))
    return energy_metadata


def get_chunk_wcs(bp, obs_id):
    """
    Set the energy WCS for the given observation.

    :param bp: The blueprint.
    :param obs_id: The current Observation ID.
    """
    logging.debug('Begin get_chunk_wcs')
    try:
        global obs_metadata
        obs_metadata = get_obs_metadata(obs_id)

        # if types contains 'AZEL_TARGET' do not create
        # spatial WCS
        # types = obs_metadata['types']
        # if 'AZEL_TARGET' not in types:
        #     bp.configure_position_axes((1, 2))

        energy_metadata = get_energy_metadata()

        # No energy metadata found
        if not energy_metadata['energy']:
            logging.error('No energy metadata found for {}'.format(obs_id))
            return

        bp.configure_energy_axis(4)
        filter_name = energy_metadata['filter_name']
        resolving_power = energy_metadata['resolving_power']
        ctype = energy_metadata['wavelength_type']
        # cunit = energy_metadata['wavelength_unit']
        naxis = energy_metadata['number_pixels']
        crpix = energy_metadata['reference_pixel']
        crval = energy_metadata['reference_wavelength']
        cdelt = energy_metadata['delta']

        # don't set the cunit since fits2caom2 sets the cunit
        # based on the ctype.
        bp.set('Chunk.energy.bandpassName', filter_name)
        bp.set('Chunk.energy.resolvingPower', resolving_power)
        bp.set('Chunk.energy.specsys', 'TOPOCENT')
        bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
        bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')
        bp.set('Chunk.energy.axis.axis.ctype', ctype)
        bp.set('Chunk.energy.axis.function.naxis', naxis)
        bp.set('Chunk.energy.axis.function.delta', cdelt)
        bp.set('Chunk.energy.axis.function.refCoord.pix', crpix)
        bp.set('Chunk.energy.axis.function.refCoord.val', crval)
    except Exception as e:
        logging.error(e)
        raise mc.CadcException(
            'Could not get chunk metadata for {}'.format(obs_id))
    logging.debug('End get_chunk_wcs')


def get_time_delta(header):
    """
    Calculate the Time WCS delta.

    :param header: The FITS header for the current extension.
    :return: The Time delta, or None if none found.
    """
    exptime = None
    if 'EXPTIME' in header:
        exptime = header['EXPTIME']
    if exptime is None:
        return None
    return float(exptime) / (24.0 * 3600.0)


def get_time_crval(header):
    """
    Calculate the Time WCS reference value.

    :param header: The FITS header for the current extension.
    :return: The Time reference value, or None if none found.
    """
    dateobs = None
    timeobs = None
    if 'DATE-OBS' in header:
        dateobs = header['DATE-OBS']
    if 'TIME-OBS' in header:
        timeobs = header['TIME-OBS']
    if not dateobs and not timeobs:
        return None
    return ac.get_datetime('{}T{}'.format(dateobs, timeobs))


def get_calibration_level(header):
    # TODO
    return CalibrationLevel.RAW_STANDARD


def get_art_product_type(header):
    obs_type = header.get('OBSTYPE')
    obs_class = header.get('OBSCLASS')

    if obs_class is None:
        global obs_metadata
        obs_class = obs_metadata['observation_class']

    logging.debug('type is {} and class is {}'.format(obs_type, obs_class))
    if obs_type is not None and obs_type == 'MASK':
        result = ProductType.AUXILIARY
    elif obs_class is None:
        if obs_type is not None and obs_type == 'OBJECT':
            obs_id = header.get('DATALAB')
            if obs_id is not None and 'CAL' in obs_id:
                result = ProductType.CALIBRATION
            else:
                result = ProductType.SCIENCE
        else:
            result = ProductType.CALIBRATION
    elif obs_class is not None and obs_class == 'science':
        result = ProductType.SCIENCE
    else:
        result = ProductType.CALIBRATION
    return result


def get_data_product_type(header):
    global obs_metadata
    mode = mc.response_lookup(obs_metadata, 'mode')
    obs_type = header.get('OBSTYPE')
    if ((mode is not None and mode == 'imaging') or
            (obs_type is not None and obs_type == 'MASK')):
        result = DataProductType.IMAGE
    else:
        result = DataProductType.SPECTRUM
    return result


def get_exposure(header):
    """
    Calculate the exposure time.

    :param header:  The FITS header for the current extension.
    :return: The exposure time, or None if not found.
    """
    if 'EXPTIME' in header:
        return header['EXPTIME']
    return None


def get_obs_intent(header):
    """
    Determine the Observation intent.

    :param header:  The FITS header for the current extension.
    :return: The Observation intent, or None if not found.
    """
    lookup = header.get('OBSCLASS')
    if lookup is None:
        object = header.get('OBJECT')
        if object in ['GCALflat', 'Bias', 'Twilight']:
            result = ObservationIntentType.CALIBRATION
        else:
            result = ObservationIntentType.SCIENCE
    elif 'science' in lookup:
        result = ObservationIntentType.SCIENCE
    else:
        result = ObservationIntentType.CALIBRATION
    return result


def get_obs_type(header):
    """
    Determine the Observation type.

    :param header:  The FITS header for the current extension.
    :return: The Observation type, or None if not found.
    """
    result = header.get('OBSTYPE')
    obs_class = header.get('OBSCLASS')
    if obs_class is not None and 'acq' in obs_class:
        result = 'ACQUISITION'
    return result


def get_target_type(header):
    global obs_metadata
    spectroscopy = mc.response_lookup(obs_metadata, 'spectroscopy')
    if spectroscopy:
        return TargetType.OBJECT
    else:
        return TargetType.FIELD


def accumulate_fits_bp(bp, uri, obs_id):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug('Begin accumulate_fits_bp.')

    bp.set('Observation.intent', 'get_obs_intent(header)')
    bp.set('Observation.type', 'get_obs_type(header)')

    bp.clear('Observation.metaRelease')
    bp.add_fits_attribute('Observation.metaRelease', 'DATE')

    bp.set('Observation.target.type', 'get_target_type(header)')

    bp.set('Plane.dataProductType', 'get_data_product_type(header)')
    bp.set('Plane.calibrationLevel', 'get_calibration_level(header)')

    bp.clear('Plane.metaRelease')
    bp.add_fits_attribute('Plane.metaRelease', 'DATE')

    bp.set('Artifact.productType', 'get_art_product_type(header)')

    bp.configure_position_axes((1, 2))
    bp.configure_time_axis(3)

    bp.set('Chunk.time.resolution', 'get_exposure(header)')
    bp.set('Chunk.time.exposure', 'get_exposure(header)')

    bp.set('Chunk.time.axis.axis.ctype', 'TIME')
    bp.set('Chunk.time.axis.axis.cunit', 'd')
    bp.set('Chunk.time.axis.error.syser', '1e-07')
    bp.set('Chunk.time.axis.error.rnder', '1e-07')
    bp.set('Chunk.time.axis.function.naxis', '1')
    bp.set('Chunk.time.axis.function.delta', 'get_time_delta(header)')
    bp.set('Chunk.time.axis.function.refCoord.pix', '0.5')
    bp.add_fits_attribute('Chunk.time.axis.function.refCoord.val', 'MJD-OBS')

    get_chunk_wcs(bp, obs_id)

    logging.debug('Done accumulate_fits_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes, must
    have this signature for import_module loading and execution.

    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.error('Begin update.')
    mc.check_param(observation, Observation)

    # for p in observation.planes:
    #     for a in observation.planes[p].artifacts:
    #         for part in observation.planes[p].artifacts[a].parts:
    #             if observation.planes[p].artifacts[a].parts[part].name == '0':
    #                 observation.planes[p].artifacts[a].parts[part].chunks.pop()
    #                 logging.error('Set chunks to None for 0-th part.')
    logging.error('Done update.')
    return True


def _build_blueprints(uris, obs_id):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param uris The list of artifact URIs for the files to be processed.
    :param obs_id The Observation ID of the file."""
    module = importlib.import_module(__name__)
    blueprints = {}
    for uri in uris:
        blueprint = ObsBlueprint(module=module)
        if not GemName.is_preview(uri):
            accumulate_fits_bp(blueprint, uri, obs_id)
        blueprints[uri] = blueprint
    return blueprints


def _get_uris(args):
    result = []
    if args.local:
        for ii in args.local:
            result.append(GemName(
                fname_on_disk=os.path.basename(ii)).file_uri)
    elif args.lineage:
        for ii in args.lineage:
            ignore, temp = mc.decompose_lineage(ii)
            result.append(temp)
    else:
        raise mc.CadcException(
            'Could not define uri from these args {}'.format(args))
    return result


def _get_obs_id(args):
    result = None
    if args.lineage:
        for lineage in args.lineage:
            temp = lineage.split('/', 1)
            if temp[1].endswith('.jpg'):
                pass
            else:
                result = temp[0]
    # TODO - figure out what to do about GemName.obs_id value
    # elif args.local:
    #     for local in args.local:
    #         if local.endswith('.jpg'):
    #             pass
    #         else:
    #             # result = GemName(
    #             #     fname_on_disk=os.path.basename(local))._get_obs_id()
    else:
        raise mc.CadcException(
            'Cannot get the obsID without the file_uri from args {}'
                .format(args))
    return result


def main_app2():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        uris = _get_uris(args)
        obs_id = _get_obs_id(args)
        blueprints = _build_blueprints(uris, obs_id)
        gen_proc(args, blueprints)
    except Exception as e:
        logging.error('Failed {} execution for {}.'.format(APPLICATION, args))
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)

    logging.debug('Done {} processing.'.format(APPLICATION))
