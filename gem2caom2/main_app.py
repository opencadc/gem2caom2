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
import re
import traceback

from caom2 import Observation, ObservationIntentType, DataProductType
from caom2 import CalibrationLevel, TargetType, ProductType, Chunk
from caom2 import SpatialWCS, CoordAxis2D, Axis, CoordPolygon2D, ValueCoord2D
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import manage_composable as mc
from caom2pipe import execute_composable as ec
from caom2pipe import astro_composable as ac

import gem2caom2.external_metadata as em
from gem2caom2.svofps import filter_metadata


__all__ = ['main_app2', 'update', 'GemName', 'COLLECTION', 'APPLICATION',
           'SCHEME', 'ARCHIVE']


APPLICATION = 'gem2caom2'
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
        # TODO - how important is the file name case? check with DB.
        # return GemName.remove_extensions(file_name.lower()).upper()
        return GemName.remove_extensions(file_name)

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        return name.replace('.fits', '').replace('.gz', ''). \
            replace('.header', '').replace('.jpg', '')

    @staticmethod
    def is_preview(entry):
        return '.jpg' in entry


def get_niri_filter_name(header):
    """
    Create the filter name for NIRI.

    :param header: The FITS header for the current extension.
    :return: The NIRI filter name, or None if none found.
    """
    filters = None
    header_filters = []
    filters2ignore = ['open', 'INVALID', 'PK50', 'pupil']
    for key in header.keys():
        if 'FILTER' in key:
            if any(x in key for x in filters2ignore):
                continue
            else:
                filtr = "".join(re.findall(r'\'(.+?)\'', header))
                filtr = filtr.replace('_', '-').strip()
                filtr = ''.join('' if ch in '()' else ch for ch in filtr)
                header_filters.append(filtr)
        filters = "&".join(header_filters)
    if filters:
        filters = re.sub(r'&', ' & ', filters)
        filters = re.sub(r'-G.{4}(|w)', '', filters)
    return filters


def get_energy_metadata(file_id):
    """
    For the given observation retrieve the energy metadata.

    :return: Dictionary of energy metadata.
    """
    logging.debug('Begin get_energy_metadata')
    instrument = em.om.get('instrument')
    if instrument in ['GMOS-N', 'GMOS-S']:
        energy_metadata = em.gmos_metadata()
    elif instrument in ['NIRI']:
        energy_metadata = em.niri_metadata(file_id)
    else:
        energy_metadata = {'energy': False}
        # raise mc.CadcException(
        #     'Do not understand energy for instrument {}'.format(instrument))
    logging.debug(
        'End get_energy_metadata for instrument {}'.format(instrument))
    return energy_metadata


def get_chunk_wcs(bp, obs_id, file_id):
    """
    Set the energy WCS for the given observation.

    :param bp: The blueprint.
    :param obs_id: The Observation ID.
    :param file_id: The file ID.
    """
    logging.debug('Begin get_chunk_wcs')
    try:

        # if types contains 'AZEL_TARGET' do not create spatial WCS
        # types = obs_metadata['types']
        # if 'AZEL_TARGET' not in types:
        #     bp.configure_position_axes((1, 2))

        energy_metadata = get_energy_metadata(file_id)

        # No energy metadata found
        if not energy_metadata['energy']:
            logging.error('No energy metadata found for '
                          'obs id {}, file id {}'.format(obs_id, file_id))
            return

        bp.configure_energy_axis(4)
        filter_name = mc.response_lookup(energy_metadata, 'filter_name')
        resolving_power = mc.response_lookup(
            energy_metadata, 'resolving_power')
        ctype = mc.response_lookup(energy_metadata, 'wavelength_type')
        naxis = mc.response_lookup(energy_metadata, 'number_pixels')
        crpix = mc.response_lookup(energy_metadata, 'reference_pixel')
        crval = mc.response_lookup(energy_metadata, 'reference_wavelength')
        cdelt = mc.response_lookup(energy_metadata, 'delta')

        # don't set the cunit since fits2caom2 sets the cunit
        # based on the ctype.
        instrument = em.om.get('instrument')
        if instrument is not None and instrument in ['NIRI']:
            bp.set('Chunk.energy.bandpassName', 'get_niri_filter_name(header)')
        else:
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
    exptime = header.get('EXPTIME')
    if exptime is None:
        return None
    return float(exptime) / (24.0 * 3600.0)


def get_time_crval(header):
    """
    Calculate the Time WCS reference value.

    :param header: The FITS header for the current extension.
    :return: The Time reference value, or None if none found.
    """
    date_obs = header.get('DATE-OBS')
    time_obs = header.get('TIME-OBS')
    if not date_obs and not time_obs:
        return None
    return ac.get_datetime('{}T{}'.format(date_obs, time_obs))


def get_calibration_level(header):
    reduction = em.om.get('reduction')
    result = CalibrationLevel.RAW_STANDARD
    if reduction is not None and 'PROCESSED_SCIENCE' in reduction:
        result = CalibrationLevel.CALIBRATED
    return result


def get_art_product_type(header):
    """
    Calculate the Artifact ProductType.

    If obsclass is unknown then CAOM2 ProductType is set to CALIBRATION.
    This should only effect early data when OBSCLASS was not in the JSON
    summary metadata or FITS headers.

    :param header:  The FITS header for the current extension.
    :return: The Artifact ProductType, or None if not found.
    """
    obs_type = _get_obs_type(header)
    obs_class = _get_obs_class(header)

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
    """
    Calculate the Plane DataProductType.

    :param header:  The FITS header for the current extension.
    :return: The Plane DataProductType, or None if not found.
    """
    mode = em.om.get('mode')
    obs_type = _get_obs_type(header)
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
    return header.get('EXPTIME')


def get_meta_release(header):
    """
    Determine the metadata release date (Observation and Plane-level).

    :param header:  The FITS header for the current extension.
    :return: The Observation/Plane release date, or None if not found.
    """
    meta_release = header.get('DATE-OBS')
    if meta_release is None:
        meta_release = em.om.get('release')
    return meta_release


def get_obs_intent(header):
    """
    Determine the Observation intent.

    :param header:  The FITS header for the current extension.
    :return: The Observation intent, or None if not found.
    """
    result = ObservationIntentType.CALIBRATION
    cal_values = ['GCALflat', 'Bias', 'BIAS', 'Twilight', 'Ar', 'FLAT', 'ARC']
    lookup = _get_obs_class(header)
    if lookup is None:
        object_value = header.get('OBJECT')
        if object_value is not None:
            instrument = header.get('INSTRUME')
            if instrument is not None and instrument == 'GRACES':
                obs_type = _get_obs_type(header)
                if obs_type is not None and obs_type in cal_values:
                    result = ObservationIntentType.CALIBRATION
                else:
                    result = ObservationIntentType.SCIENCE
            else:
                if object_value in cal_values:
                    result = ObservationIntentType.CALIBRATION
                else:
                    result = ObservationIntentType.SCIENCE
    elif 'science' in lookup:
        result = ObservationIntentType.SCIENCE
    return result


def get_obs_type(header):
    """
    Determine the Observation type.

    :param header:  The FITS header for the current extension.
    :return: The Observation type, or None if not found.
    """
    result = _get_obs_type(header)
    obs_class = _get_obs_class(header)
    if obs_class is not None and (obs_class == 'acq' or obs_class == 'acqCal'):
        result = 'ACQUISITION'
    return result


def get_target_type(header):
    """
    Calculate the Target TargetType

    :param header:  The FITS header for the current extension.
    :return: The Target TargetType, or None if not found.
    """
    spectroscopy = em.om.get('spectroscopy')
    if spectroscopy:
        return TargetType.OBJECT
    else:
        return TargetType.FIELD


def _get_obs_class(header):
    """Common location to lookup OBSCLASS from the FITS headers, and if
    it's not present, to lookup observation_class from JSON summary
    metadata."""
    obs_class = header.get('OBSCLASS')
    if obs_class is None:
        obs_class = em.om.get('observation_class')
    return obs_class


def _get_obs_type(header):
    """Common location to lookup OBSTYPE from the FITS headers, and if
    it's not present, to lookup observation_type from JSON summary
    metadata."""
    obs_type = header.get('OBSTYPE')
    if obs_type is None:
        obs_type = em.om.get('observation_type')
    return obs_type


def accumulate_fits_bp(bp, obs_id, file_id):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug('Begin accumulate_fits_bp.')
    em.get_obs_metadata(file_id)

    bp.set('Observation.intent', 'get_obs_intent(header)')
    bp.set('Observation.type', 'get_obs_type(header)')
    bp.set('Observation.metaRelease', 'get_meta_release(header)')
    bp.set('Observation.target.type', 'get_target_type(header)')
    # GRACES has entries with RUNID set to non-proposal ID values
    # so clear the default FITS value lookup
    bp.clear('Observation.proposal.id')

    bp.set('Plane.dataProductType', 'get_data_product_type(header)')
    bp.set('Plane.calibrationLevel', 'get_calibration_level(header)')
    bp.set('Plane.metaRelease', 'get_meta_release(header)')

    bp.set('Plane.provenance.name', 'Gemini Observatory Data')
    bp.set('Plane.provenance.project', 'Gemini Archive')
    bp.add_fits_attribute('Plane.provenance.producer', 'IMAGESWV')
    bp.set_default('Plane.provenance.producer', 'Gemini Observatory')
    bp.set('Plane.provenance.reference',
           'http://archive.gemini.edu/searchform/{}'.format(obs_id))

    bp.set('Artifact.productType', 'get_art_product_type(header)')

    bp.configure_position_axes((1, 2))
    bp.configure_time_axis(3)

    # The Chunk time metadata is calculated using keywords from the
    # primary header, and the only I could figure out to access keywords
    # in the primary is through a function. JB
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

    get_chunk_wcs(bp, obs_id, file_id)

    logging.debug('Done accumulate_fits_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes, must
    have this signature for import_module loading and execution.

    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.error('Begin update.')
    mc.check_param(observation, Observation)

    headers = None
    if 'headers' in kwargs:
        headers = kwargs['headers']

    for p in observation.planes:
        mode = observation.planes[p].data_product_type
        for a in observation.planes[p].artifacts:
            for part in observation.planes[p].artifacts[a].parts:
                for c in observation.planes[p].artifacts[a].parts[part].chunks:
                    if observation.instrument.name == 'NIRI':
                        _update_chunk_energy(c, mode, headers)

                    if (observation.instrument.name == 'GRACES' and
                            mode == DataProductType.SPECTRUM):
                        _update_chunk_position(c)
    logging.error('Done update.')
    return True


def _update_chunk_energy(chunk, mode, headers):
    """Check that position information has been set appropriately for
    spectra. Guidance given is it must be a bounds."""
    logging.debug('Begin _update_chunk_position')
    mc.check_param(chunk, Chunk)

    # No energy information is determined for darks.  The
    # latter are sometimes only identified by a 'blank' filter.  e.g.
    # NIRI 'flats' are sometimes obtained with the filter wheel blocked off.

    bandpass_name = chunk.energy.bandpass_name
    if 'blank' in bandpass_name:
        chunk.energy = None
    else:
        cval = 0.0
        delta = 0.0
        resolving_power = 0.0
        filter_md = filter_metadata(
            em.om.get('instrument'), bandpass_name)
        if mode == 'imaging':
            delta = filter_md['wl_eff_width']
            cval = filter_md['wl_eff']
            resolving_power = cval/delta
            cval /= 1.0e10
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

        chunk.energy.resolvingPower = resolving_power
        chunk.energy.specsys = 'TOPOCENT'
        chunk.energy.ssysobs = 'TOPOCENT'
        chunk.energy.ssyssrc = 'TOPOCENT'
        chunk.energy.axis.axis.ctype = 'WAVE'
        chunk.energy.axis.function.naxis = 1024
        chunk.energy.axis.function.delta = delta
        chunk.energy.axis.function.refCoord.pix = 512.0
        chunk.energy.axis.function.refCoord.val = cval


def _update_chunk_position(chunk):
    """Check that position information has been set appropriately for
    GRACES spectra. Guidance given from DB, Jan 18/19 is it must be a
    bounds."""

    logging.debug('Begin _update_chunk_position')
    mc.check_param(chunk, Chunk)

    ra = em.om.get('ra')
    dec = em.om.get('dec')

    if ra is not None and dec is not None:
        axis1 = Axis(ctype='RA---TAN', cunit=None)
        axis2 = Axis(ctype='DEC--TAN', cunit=None)
        polygon = CoordPolygon2D()
        # radius = 5.0 degrees, said DB in test data list of Jan 18/19
        for x, y in ([0, 1], [1, 1], [1, 0], [0, 0]):
            ra_pt = ra - 5.0*(0.5-float(x))
            dec_pt = dec - 5.0*(0.5-float(y))
            polygon.vertices.append(ValueCoord2D(ra_pt, dec_pt))
        polygon.vertices.append(ValueCoord2D(ra, dec))
        axis = CoordAxis2D(axis1=axis1, axis2=axis2,
                           error1=None, error2=None,
                           range=None, bounds=polygon, function=None)
        chunk.position = SpatialWCS(axis)

    logging.debug('End _update_chunk_position')


def _build_blueprints(uris, obs_id, file_id):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param uris The list of artifact URIs for the files to be processed.
    :param obs_id The Observation ID of the file.
    :param file_id The file ID."""
    module = importlib.import_module(__name__)
    blueprints = {}
    for uri in uris:
        blueprint = ObsBlueprint(module=module)
        if not GemName.is_preview(uri):
            logging.info('Building blueprint for {}'.format(uri))
            accumulate_fits_bp(blueprint, obs_id, file_id)
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
    if args.in_obs_xml:
        temp = args.in_obs_xml.name.split('/')
        result = temp[-1].split('.')[0]
    elif args.lineage:
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


def _get_file_id(args):
    result = None
    if args.lineage:
        for lineage in args.lineage:
            temp = lineage.split(':', 1)
            if temp[1].endswith('.jpg'):
                pass
            else:
                result = re.split(r"[/.]", temp[1])[1]
    else:
        raise mc.CadcException(
            'Cannot get the fileID without the file_uri from args {}'
                .format(args))
    return result


def main_app2():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        uris = _get_uris(args)
        obs_id = _get_obs_id(args)
        file_id = _get_file_id(args)
        blueprints = _build_blueprints(uris, obs_id, file_id)
        gen_proc(args, blueprints)
    except Exception as e:
        logging.error('Failed {} execution for {}.'.format(APPLICATION, args))
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)

    logging.debug('Done {} processing.'.format(APPLICATION))
