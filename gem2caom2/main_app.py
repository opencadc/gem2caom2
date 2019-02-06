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
from caom2 import SpectralWCS, CoordAxis1D, CoordFunction1D, RefCoord
from caom2 import TypedList
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import manage_composable as mc
from caom2pipe import astro_composable as ac

import gem2caom2.external_metadata as em
from gem2caom2.gem_name import GemName, COLLECTION, ARCHIVE, SCHEME

__all__ = ['main_app2', 'update', 'COLLECTION', 'APPLICATION', 'SCHEME',
           'ARCHIVE']


APPLICATION = 'gem2caom2'


def get_energy_metadata(file_id):
    """
    For the given observation retrieve the energy metadata.

    :return: Dictionary of energy metadata.
    """
    logging.debug('Begin get_energy_metadata')
    instrument = em.om.get('instrument')
    if instrument in ['GMOS-N', 'GMOS-S']:
        energy_metadata = em.gmos_metadata()
    else:
        energy_metadata = {'energy': False}
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

        _get_chunk_energy(bp, obs_id, file_id)
    except Exception as e:
        logging.error(e)
        raise mc.CadcException(
            'Could not get chunk metadata for {}'.format(obs_id))
    logging.debug('End get_chunk_wcs')


def _get_chunk_energy(bp, obs_id, file_id):
    energy_metadata = get_energy_metadata(file_id)

    # No energy metadata found
    if energy_metadata['energy']:
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
    # else:
    #     logging.info('No energy metadata found for '
    #                  'obs id {}, file id {}'.format(obs_id, file_id))


def get_time_delta(header):
    """
    Calculate the Time WCS delta.

    :param header: The FITS header for the current extension.
    :return: The Time delta, or None if none found.
    """
    exptime = get_exposure(header)
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

    # logging.error('type is {} and class is {}'.format(obs_type, obs_class))
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
            instrument = header.get('INSTRUME')
            if instrument is not None:
                if instrument == 'PHOENIX':
                    if _is_phoenix_calibration(header):
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                elif instrument == 'OSCIR':
                    if _is_oscir_calibration(header):
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                elif instrument == 'FLAMINGOS':
                    if _is_flamingos_calibration(header):
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                elif instrument == 'Hokupaa+QUIRC':
                    if _is_hokupaa_calibration(header):
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                else:
                    result = ProductType.CALIBRATION
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


def get_data_release(header):
    """
    Determine the plane-level data release date.

    :param header:  The FITS header for the current extension.
    :return: The Plane release date, or None if not found.
    """
    # every instrument has a 'release' keyword in the JSON summary
    # not every instrument (Michelle) has a RELEASE keyword in
    # the appropriate headers
    return em.om.get('release')


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
    cal_values = ['GCALflat', 'Bias', 'BIAS', 'Twilight', 'Ar', 'FLAT', 'ARC',
                  'Domeflat']
    lookup = _get_obs_class(header)
    # logging.error('lookup is {}'.format(lookup))
    if lookup is None:
        object_value = header.get('OBJECT')
        # logging.error('object_value is {}'.format(object_value))
        if object_value is not None:
            instrument = header.get('INSTRUME')
            # logging.error('instrument is {}'.format(instrument))
            if (instrument is not None and
                    instrument in ['GRACES', 'TReCS', 'PHOENIX', 'OSCIR',
                                   'Hokupaa+QUIRC']):
                if instrument == 'GRACES':
                    obs_type = _get_obs_type(header)
                    if obs_type is not None and obs_type in cal_values:
                        result = ObservationIntentType.CALIBRATION
                    else:
                        result = ObservationIntentType.SCIENCE
                elif instrument == 'TReCS':
                    data_label = header.get('DATALAB')
                    if data_label is not None and '-CAL' in data_label:
                        result = ObservationIntentType.CALIBRATION
                elif instrument == 'PHOENIX':
                    if _is_phoenix_calibration(header):
                        result = ObservationIntentType.CALIBRATION
                    else:
                        result = ObservationIntentType.SCIENCE
                elif instrument == 'OSCIR':
                    if _is_oscir_calibration(header):
                        result = ObservationIntentType.CALIBRATION
                    else:
                        result = ObservationIntentType.SCIENCE
                elif instrument == 'Hokupaa+QUIRC':
                    if _is_hokupaa_calibration(header):
                        result = ObservationIntentType.CALIBRATION
                    else:
                        result = ObservationIntentType.SCIENCE
                else:
                    logging.error(
                        'get_obs_intent: no handling for {}'.format(
                            instrument))
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


def get_target_moving(header):
    """
    Calculate whether the Target moving.
    Non-sidereal tracking -> setting moving target to "True"

    :param header:  The FITS header for the current extension.
    :return: The Target TargetType, or None if not found.
    """
    types = em.om.get('types')
    if 'NON_SIDEREAL' in types:
        return True
    else:
        return None


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
    obs_type = em.om.get('observation_type')
    if obs_type is None:
        obs_type = header.get('OBSTYPE')
    return obs_type


def _is_flamingos_calibration(header):
    object_value = em.om.get('object')
    if ('Dark' in object_value or
        'DARK' in object_value or
            'Flat' in object_value):
        return True
    return False


def _is_hokupaa_calibration(header):
    image_type = header.get('IMAGETYP')
    # took the list from the AS GEMINI Obs. Type values
    return image_type in ['DARK', 'BIAS', 'CAL', 'ARC', 'FLAT']


def _is_oscir_calibration(header):
    # TODO - don't know what else to do right now
    return False


def _is_phoenix_calibration(header):
    object_value = header.get('OBJECT').lower()
    # logging.error('object_value is {}'.format(object_value))
    if ('flat ' in object_value or
        'dark ' in object_value or
        'arc' in object_value or
        'comp' in object_value or
        'lamp' in object_value or
            'comparison' in object_value):
        return True
    return False


def accumulate_fits_bp(bp, obs_id, file_id):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug('Begin accumulate_fits_bp.')
    em.get_obs_metadata(file_id)

    bp.set('Observation.intent', 'get_obs_intent(header)')
    bp.set('Observation.type', 'get_obs_type(header)')
    bp.set('Observation.metaRelease', 'get_meta_release(header)')
    bp.set('Observation.target.type', 'get_target_type(header)')
    bp.set('Observation.target.moving', 'get_target_moving(header)')
    # GRACES has entries with RUNID set to non-proposal ID values
    # so clear the default FITS value lookup
    bp.clear('Observation.proposal.id')

    bp.set('Plane.dataProductType', 'get_data_product_type(header)')
    bp.set('Plane.calibrationLevel', 'get_calibration_level(header)')
    bp.set('Plane.metaRelease', 'get_meta_release(header)')
    bp.set('Plane.dataRelease', 'get_data_release(header)')

    bp.set('Plane.provenance.name', 'Gemini Observatory Data')
    bp.set('Plane.provenance.project', 'Gemini Archive')
    # Add IMAGESWV for GRACES
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

    try:
        for p in observation.planes:
            mode = observation.planes[p].data_product_type
            for a in observation.planes[p].artifacts:
                for part in observation.planes[p].artifacts[a].parts:

                    if part == '2' and observation.instrument.name == 'GPI':
                        # GPI datasets have two extensions. First is science
                        # image (with WCS), second is data quality for each
                        # pixel (no WCS).
                        logging.info(
                            'GPI: Setting chunks to None for part {}'.format(
                                part))
                        observation.planes[p].artifacts[a].parts[part].chunks \
                            = TypedList(Chunk,)
                        continue
                    for c in observation.planes[p].artifacts[a].parts[
                            part].chunks:

                        # energy WCS
                        if _reset_energy(headers[0].get('DATALAB')):
                            c.energy = None
                            c.energy_axis = None
                        else:
                            if observation.instrument.name == 'NIRI':
                                _update_chunk_energy_niri(c, headers)
                            elif observation.instrument.name == 'GPI':
                                _update_chunk_energy_gpi(c, headers[0])
                            elif observation.instrument.name == 'F2':
                                _update_chunk_energy_f2(c, headers)
                            elif observation.instrument.name == 'GSAOI':
                                _update_chunk_energy_gsaoi(c)
                            elif observation.instrument.name == 'NICI':
                                _update_chunk_energy_nici(c)

                        # position WCS
                        if (observation.instrument.name == 'GRACES' and
                                mode == DataProductType.SPECTRUM):
                            _update_chunk_position(c)
    except Exception as e:
        logging.error(e)
    logging.error('Done update.')
    return True


def _update_chunk_energy_niri(chunk, headers):
    """NIRI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_niri')
    mc.check_param(chunk, Chunk)

    # No energy information is determined for darks.  The
    # latter are sometimes only identified by a 'blank' filter.  e.g.
    # NIRI 'flats' are sometimes obtained with the filter wheel blocked off.

    n_axis = 1
    try:
        header = headers[0]

        filters = get_filter_name(header)
        filter_md = em.get_filter_metadata('NIRI', filters)
        filter_name = em.om.get('filter_name')

        mode = em.om.get('mode')
        if mode == 'imaging':
            logging.debug('SpectralWCS: NIRI imaging mode.')
            c_val, delta, resolving_power = _imaging_energy(filter_md)
        elif mode in ['LS', 'spectroscopy']:
            logging.debug('SpectralWCS: NIRI LS|Spectroscopy mode.')
            c_val = 0.0  # TODO
            delta = filter_md['wl_eff_width'] / n_axis / 1.0e10
            fp_mask = header.get('FPMASK')
            bandpass_name = filter_name[0]
            f_ratio = fp_mask.split('_')[0]
            logging.debug('Bandpass name is {} f_ratio is {}'.format(
                bandpass_name, f_ratio))
            if bandpass_name in em.NIRI_RESOLVING_POWER:
                resolving_power = \
                    em.NIRI_RESOLVING_POWER[bandpass_name][f_ratio]
            else:
                resolving_power = None
                logging.debug('No resolving power.')
        else:
            raise mc.CadcException(
                'Do not understand mode {}'.format(mode))

        _build_chunk_energy(chunk, n_axis, c_val, delta, filter_name,
                            resolving_power)
    except:
        tb = traceback.format_exc()
        logging.error(tb)
    logging.debug('End _update_chunk_energy_niri')


def _build_chunk_energy(chunk, n_axis, c_val, delta, filter_name,
                        resolving_power):
    # build the CAOM2 structure

    # DB - 02-04-19 - initial pass, do not try to calculate
    # dispersion, so naxis=1, prefer um as units, strip the bar code
    # from the filter names
    chunk.energy_axis = 4
    axis = Axis(ctype='WAVE', cunit='um')
    ref_coord = RefCoord(pix=float(n_axis/2.0), val=c_val)
    # SGo - assume a function until DB says otherwise
    function = CoordFunction1D(naxis=n_axis,
                               delta=delta,
                               ref_coord=ref_coord)
    coord_axis = CoordAxis1D(axis=axis,
                             error=None,
                             range=None,
                             bounds=None,
                             function=function)
    chunk.energy = SpectralWCS(axis=coord_axis,
                               specsys='TOPOCENT',
                               ssysobs='TOPOCENT',
                               ssyssrc='TOPOCENT',
                               restfrq=None,
                               restwav=None,
                               velosys=None,
                               zsource=None,
                               velang=None,
                               bandpass_name=filter_name,
                               transition=None,
                               resolving_power=resolving_power)


def _update_chunk_energy_gpi(chunk, header):
    """NIRI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_gpi')
    mc.check_param(chunk, Chunk)

    n_axis = 1

    try:
        filter_name = em.om.get('filter_name')
        filter_md = em.get_filter_metadata('GPI', filter_name)

        mode = em.om.get('mode')
        if mode in ['imaging', 'IFP', 'IFS']:
            logging.debug('SpectralWCS: GPI imaging mode.')
            c_val, delta, resolving_power = _imaging_energy(filter_md)
        elif mode in ['LS', 'spectroscopy']:
            logging.debug('SpectralWCS: GPI LS|Spectroscopy mode.')
            c_val = 0.0  # TODO
            delta = filter_md['wl_eff_width'] / n_axis / 1.0e10
            fp_mask = header.get('FPMASK')
            bandpass_name = filter_name[0]
            f_ratio = fp_mask.split('_')[0]
            logging.debug('Bandpass name is {} f_ratio is {}'.format(
                bandpass_name, f_ratio))
            if bandpass_name in em.NIRI_RESOLVING_POWER:
                resolving_power = \
                    em.NIRI_RESOLVING_POWER[bandpass_name][f_ratio]
            else:
                resolving_power = None
                logging.debug('No resolving power.')
        else:
            raise mc.CadcException(
                'Do not understand mode {}'.format(mode))

        _build_chunk_energy(chunk, n_axis, c_val, delta, filter_name,
                            resolving_power)
    except:
        tb = traceback.format_exc()
        logging.error(tb)
    logging.debug('End _update_chunk_energy_gpi')


def _update_chunk_energy_f2(chunk, headers):
    """NIRI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_f2')
    mc.check_param(chunk, Chunk)

    n_axis = 1
    try:
        header = headers[0]

        filter_name = em.om.get('filter_name')
        filter_md = em.get_filter_metadata('Flamingos2', filter_name)

        mode = em.om.get('mode')
        logging.error('mode is {}'.format(mode))
        if mode in ['imaging', 'IFP', 'IFS']:
            logging.debug('SpectralWCS: F2 imaging mode.')
            reference_wavelength, delta, resolving_power = \
                _imaging_energy(filter_md)
        elif mode in ['LS', 'spectroscopy', 'MOS']:
            logging.debug('SpectralWCS: F2 LS|Spectroscopy mode.')
            fp_mask = header.get('MASKNAME')
            if mode == 'LS':
                slit_width = fp_mask[0]
            else:
                # DB - 04-04-19 no way to determine slit widths used in
                # custom mask, so assume 2
                slit_width = '2'
            n_axis = 2048
            delta = filter_md['wl_eff_width'] / n_axis / 1.0e4
            reference_wavelength = em.om.get('central_wavelength')
            grism_name = header.get('GRISM')
            logging.error('grism name is {} fp_mask is {}'.format(grism_name, fp_mask))
            # lookup values from
            # https://www.gemini.edu/sciops/instruments/flamingos2/spectroscopy/longslit-spectroscopy
            lookup = {'1': [1300.0, 3600.0],
                      '2': [900.0, 2800.0],
                      '3': [600.0, 1600.0],
                      '4': [350.0, 1300.0],
                      '6': [130.0, 1000.0],
                      '8': [100.0, 750.0]}
            if grism_name.startswith('R3K_'):
                resolving_power = lookup[slit_width][1]
            else:
                resolving_power = lookup[slit_width][0]
        else:
            raise mc.CadcException(
                'Do not understand mode {}'.format(mode))

        _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                            filter_name, resolving_power)
    except:
        tb = traceback.format_exc()
        logging.error(tb)
    logging.debug('End _update_chunk_energy_f2')


def _update_chunk_energy_gsaoi(chunk):
    """NIRI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_gsaoi')
    mc.check_param(chunk, Chunk)

    n_axis = 1
    filter_name = em.om.get('filter_name')
    filter_md = em.get_filter_metadata('GSAOI', filter_name)

    mode = em.om.get('mode')
    logging.error('mode is {}'.format(mode))
    if mode == 'imaging':
        logging.debug('SpectralWCS: GSAOI imaging mode.')
        reference_wavelength, delta, resolving_power = \
            _imaging_energy(filter_md)
    else:
        raise mc.CadcException(
            'Do not understand mode {}'.format(mode))

    _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                        filter_name, resolving_power)

    logging.debug('End _update_chunk_energy_gsaoi')


def _update_chunk_energy_nici(chunk):
    """NICI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_nici')
    mc.check_param(chunk, Chunk)

    n_axis = 1
    filter_name = em.om.get('filter_name')
    filter_md = em.get_filter_metadata('NICI', filter_name)

    mode = em.om.get('mode')
    logging.error('mode is {}'.format(mode))
    if mode == 'imaging':
        logging.debug('SpectralWCS: NICI imaging mode.')
        reference_wavelength, delta, resolving_power = \
            _imaging_energy(filter_md)
    else:
        raise mc.CadcException(
            'Do not understand mode {}'.format(mode))

    _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                        filter_name, resolving_power)

    logging.debug('End _update_chunk_energy_nici')


def _reset_energy(data_label):
    """
    Return True if there should be no energy WCS information created at
    the chunk level.

    :param data_label str for useful logging information only.
    """
    result = False
    observation_type = em.om.get('observation_type')
    om_filter_name = em.om.get('filter_name')

    if observation_type == 'DARK' or 'blank' in om_filter_name:
        logging.info('No chunk energy for {}'.format(data_label))
        result = True
    return result


def _imaging_energy(filter_md):
    """"""
    # all the filter units are Angstroms, and the chunk-level energy WCS
    # unit of choice is 'um', therefore 1.0e4
    c_val = filter_md['wl_eff'] / 1.0e4
    delta = filter_md['wl_eff_width'] / 1.0e4
    resolving_power = c_val / delta
    return c_val, delta, resolving_power


def get_filter_name(header):
    """
    Create the filter names.

    :param header: The FITS header for the current extension.
    :return: The filter names, or None if none found.
    """
    filters = None
    header_filters = []

    # DB - 04-02-19 - strip out anything with 'pupil' as it doesn't affect
    # energy transmission
    filters2ignore = ['open', 'invalid', 'pupil']
    for key in header.keys():
        if 'FILTER' in key:
            value = header.get(key).lower()
            ignore = False
            for ii in filters2ignore:
                if ii.startswith(value):
                    ignore = True
                    break
            if ignore:
                continue
            else:
                header_filters.append(header.get(key).strip())
        filters = '&'.join(header_filters)
    logging.info('Filters are {}'.format(filters))
    return filters


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
