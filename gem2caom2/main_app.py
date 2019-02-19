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
from caom2utils import WcsParser
from caom2pipe import manage_composable as mc
from caom2pipe import astro_composable as ac

import gem2caom2.external_metadata as em
from gem2caom2.gem_name import GemName, COLLECTION, ARCHIVE, SCHEME

__all__ = ['main_app2', 'update', 'COLLECTION', 'APPLICATION', 'SCHEME',
           'ARCHIVE']

# TODO LIST:
#
# PI information
# fix phoenix spatial
# spectroscopy check
# flamingos
# hrwfs
# hokupaa
# oscir
# bHROS

APPLICATION = 'gem2caom2'
# GPI radius == 2.8 arcseconds, according to Christian Marois via DB 02-07-19
# DB - 18-02-19 - Replace “5.0” for “2.8" for GPI field of view.
# NOTE:  To be more accurate for GRACES this size could be reduced
# from 5" to 1.2" since that’s the size of the fibre.
RADIUS_LOOKUP = {'GPI': 2.8 / 3600.0,  # units are arcseconds
                 'GRACES': 1.2 / 3600.0,
                 'PHOENIX': 5.0 / 3600.0}
phoenix_spectral_axis = 1024


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


def get_cd11(header):
    cdelt1 = header.get('CDELT1')
    if cdelt1 is None:
        instrument = em.om.get('instrument')
        result = RADIUS_LOOKUP[instrument]
    else:
        result = cdelt1
    return result


def get_cd22(header):
    cdelt2 = header.get('CDELT2')
    if cdelt2 is None:
        instrument = em.om.get('instrument')
        result = RADIUS_LOOKUP[instrument]
    else:
        result = cdelt2
    return result


def get_data_product_type(header):
    """
    Calculate the Plane DataProductType.

    :param header:  The FITS header for the current extension.
    :return: The Plane DataProductType, or None if not found.
    """
    instrument = em.om.get('instrument')
    if instrument == 'FLAMINGOS':
        result, ignore = _get_flamingos_mode(header)
    else:
        mode = em.om.get('mode')
        obs_type = _get_obs_type(header)
        if mode is None:
            raise mc.CadcException('No mode information found for {}'.format(
                em.om.get('filename')))
        elif ((mode is not None and mode == 'imaging') or
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
    Calculate the exposure time. EXPTIME in the header is not always the
    total exposure time for some of the IR instruments.  For these EXPTIME
    is usually the individual exposure time for each chop/nod and a bunch of
    chops/nods are executed for one observation.  Gemini correctly
    allows for this in the json ‘exposure_time’ value.

    :param header:  The FITS header for the current extension (unused).
    :return: The exposure time, or None if not found.
    """
    return em.om.get('exposure_time')


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


def get_naxis2_phoenix(header):
    """Over-ride the NAXIS2 length, but preserve the original value for
    use later by Spectral WCS computations. This only works because the
    blueprint applies the functions before it sets the over-ride values in
    the headers."""
    global phoenix_spectral_axis
    phoenix_spectral_axis = header.get('NAXIS2')
    return 1


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


def get_time_function_val(header):
    """
    Calculate the Chunk Time WCS function value, in 'mjd'.

    :param header:  The FITS header for the current extension (not used).
    :return: The Time WCS value from JSON Summary Metadata.
    """
    instrument = em.om.get('instrument')
    if instrument == 'FLAMINGOS':
        # Another FLAMINGOS correction needed:  DATE-OBS in header and
        # json doesn’t include the time  but sets it to 00:00:00.  You have
        # two choices:  concatenate DATE-OBS and UTC header values to the
        # standard form “2002-11-06T07:06:00.3” and use as you use json
        # value for computing temporal WCS data or, for FLAMINGOS only,
        # use the MJD header value in the header as the CRVAL for time.
        time_val = header.get('MJD')
    else:
        time_string = em.om.get('ut_datetime')
        time_val = ac.get_datetime(time_string)
    return time_val


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


def _get_flamingos_mode(header):
    # DB - 18-02-19
    # Determine mode from DECKER and GRISM keywords:
    #
    # If DECKER = imaging
    # observing mode = imaging
    # else if DECKER = mos or slit
    # if GRISM contains the string ‘open’
    # observation mode = imaging and observation type = ACQUISITION
    # else if GRISM contains the string ‘dark’
    # observation mode = spectroscopy and observation type = DARK
    # else
    # observation mode = spectroscopy
    # else
    # error? DECKER might also be ‘dark’ at times in which case
    # observation type = dark and observation mode can be either
    # spectroscopy or imaging but I haven’t found an example yet.
    #
    # As indicated above in if/else code, if GRISM  keyword
    # contains substring ‘dark’ then that is another way to
    # identify a dark observation.
    decker = header.get('DECKER')
    data_type = DataProductType.SPECTRUM
    obs_type = None
    if decker is not None:
        if decker == 'imaging':
            data_type = DataProductType.IMAGE
        elif decker in ['mos', 'slit']:
            grism = header.get('GRISM')
            if grism is not None:
                if 'open' in grism:
                    data_type = DataProductType.IMAGE
                    obs_type = 'ACQUISITION'
                elif 'dark' in grism:
                    data_type = DataProductType.SPECTRUM
                    obs_type = 'DARK'
    else:
        raise mc.CadcException('No mode information found for {}'.format(
            em.om.get('filename')))
    if obs_type is None:
        # DB - Also, since I’ve found FLAMINGOS spectra if OBJECT keyword or
        # json value contains ‘arc’ (any case) then it’s an ARC observation
        # type
        #
        # For FLAMINGOS since OBS_TYPE seems to always be st to ‘Object’
        # could in principal look at the OBJECT keyword value or json value:
        # if it contains ‘flat’ as a substring (any case) then set
        # observation type to ‘flat’.  Ditto for ‘dark’
        object_value = em.om.get('object')
        if object_value is None:
            object_value = header.get('OBJECT')
        object_value = object_value.lower()
        if 'arc' in object_value:
            obs_type = 'ARC'
        elif 'flat' in object_value:
            obs_type = 'FLAT'
        elif 'dark' in object_value:
            obs_type = 'DARK'
    return data_type, obs_type


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

    instrument = em.om.get('instrument')
    mode = em.om.get('mode')
    if (instrument in ['GPI', 'PHOENIX'] or
            (instrument == 'GRACES' and mode is not None and
                mode != 'imaging')):
        # DB - 18-02-19 - for hard-coded field of views use:
        # CRVAL1  = RA value from json or header (degrees
        # CRVAL2  = Dec value from json or header (degrees)
        # CDELT1  = 5.0/3600.0 (Plate scale along axis1 in degrees/pixel
        #           for 5" size)
        # CDELT2  = 5.0/3600.0
        # CROTA1  = 0.0 / Rotation in degrees
        # NAXIS1 = 1
        # NAXIS2 = 1
        # CRPIX1 = 1.0
        # CRPIX2 = 1.0
        # CTYPE1 = RA---TAN
        # CTYPE2 = DEC--TAN
        #
        # May be slightly different for Phoenix if headers give the width
        # and rotation of the slit.  i.e CDELT1 and 2 may be different and
        # CROTA1 non-zero.

        bp.set('Chunk.position.axis.axis1.ctype', 'RA---TAN')
        bp.set('Chunk.position.axis.axis2.ctype', 'DEC--TAN')
        bp.set('Chunk.position.axis.axis1.cunit', 'deg')
        bp.set('Chunk.position.axis.axis2.cunit', 'deg')
        bp.set('Chunk.position.axis.function.dimension.naxis1', '1')
        bp.set('Chunk.position.axis.function.dimension.naxis2', '1')
        bp.set('Chunk.position.axis.function.refCoord.coord1.pix', '1.0')
        bp.set('Chunk.position.axis.function.refCoord.coord2.pix', '1.0')
        bp.clear('Chunk.position.axis.function.refCoord.coord1.val')
        bp.add_fits_attribute(
            'Chunk.position.axis.function.refCoord.coord1.val', 'RA')
        bp.clear('Chunk.position.axis.function.refCoord.coord2.val')
        bp.add_fits_attribute(
            'Chunk.position.axis.function.refCoord.coord2.val', 'DEC')
        bp.set('Chunk.position.axis.function.cd11', 'get_cd11(header)')
        bp.set('Chunk.position.axis.function.cd22', 'get_cd22(header)')
        bp.set_default('Chunk.position.axis.function.cd12', '0.0')
        bp.set_default('Chunk.position.axis.function.cd21', '0.0')

        if instrument == 'PHOENIX':
            bp.set('Chunk.position.axis.function.dimension.naxis2',
                   'get_naxis2_phoenix(header)')

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
    bp.set('Chunk.time.axis.function.refCoord.val',
           'get_time_function_val(header)')

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

                if (observation.instrument.name in
                        ['michelle', 'TReCS', 'NIFS', 'GNIRS']):
                    # Michelle is a retired visitor instrument.
                    # Spatial WCS info is in primary header. There
                    # are a variable number of FITS extensions
                    # defined by primary keyword NUMEXT; assume the
                    # same spatial WCS for each for now - it differs
                    # only slightly because of telescope 'chopping'
                    # and 'nodding' used in acquisition. DB - 01-18-19
                    #
                    # DB - 01-18-19 - NIFS has no WCS info in extension; use
                    # primary header, GNIRS has no WCS info in extension; use
                    # primary header
                    _update_position_from_zeroth_header(
                        observation.planes[p].artifacts[a], headers,
                        observation.observation_id)

                for part in observation.planes[p].artifacts[a].parts:

                    if part == '2' and observation.instrument.name == 'GPI':
                        # GPI datasets have two extensions. First is science
                        # image (with WCS), second is data quality for each
                        # pixel (no WCS).
                        logging.info(
                            'GPI: Setting chunks to None for part {} for {}'.format(
                                part, observation.observation_id))
                        observation.planes[p].artifacts[a].parts[part].chunks \
                            = TypedList(Chunk,)
                        continue
                    for c in observation.planes[p].artifacts[a].parts[
                            part].chunks:
                        header = headers[int(part)]

                        # energy WCS
                        if _reset_energy(headers[0].get('DATALAB')):
                            c.energy = None
                            c.energy_axis = None
                        else:
                            if observation.instrument.name == 'NIRI':
                                _update_chunk_energy_niri(
                                    c, headers[0], observation.observation_id)
                            elif observation.instrument.name == 'GPI':
                                _update_chunk_energy_gpi(c, headers[0])
                            elif observation.instrument.name == 'F2':
                                _update_chunk_energy_f2(c, headers[0])
                            elif observation.instrument.name == 'GSAOI':
                                _update_chunk_energy_gsaoi(c)
                            elif observation.instrument.name == 'NICI':
                                _update_chunk_energy_nici(
                                    c, header, observation.observation_id)
                            elif observation.instrument.name == 'TReCS':
                                _update_chunk_energy_trecs(
                                    c, header, observation.observation_id)
                            elif observation.instrument.name == 'michelle':
                                _update_chunk_energy_michelle(
                                    c, observation.observation_id)
                            elif observation.instrument.name == 'NIFS':
                                _update_chunk_energy_nifs(
                                    c, header, observation.observation_id)
                            elif observation.instrument.name == 'GNIRS':
                                _update_chunk_energy_gnirs(c, header)
                            elif observation.instrument.name == 'PHOENIX':
                                _update_chunk_energy_phoenix(
                                    c, header, observation.observation_id)
                            elif observation.instrument.name == 'FLAMINGOS':
                                _update_chunk_energy_flamingos(
                                    c, header,
                                    observation.planes[p].data_product_type,
                                    observation.observation_id)
                                ignore, obs_type = _get_flamingos_mode(header)
                                if (obs_type is not None and
                                        observation.type is None):
                                    observation.type = obs_type
                            # elif observation.instrument.name == 'hrwfs':
                            #     _update_chunk_energy_hrwfs(
                            #         c, header, observation.observation_id)

                        # position WCS
                        if (part == '1' and
                                observation.instrument.name == 'GPI'):
                            # equinox information only available from the
                            # 0th header
                            c.position.equinox = headers[0].get('TRKEQUIN')
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.error(tb)
        return None
    logging.error('Done update.')
    return observation


def _build_chunk_energy(chunk, n_axis, c_val, delta, filter_name,
                        resolving_power):
    # build the CAOM2 structure

    # DB - 02-04-19 - initial pass, do not try to calculate
    # dispersion, so prefer um as units, strip the bar code
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
    bandpass_name = None if len(filter_name) == 0 else filter_name
    chunk.energy = SpectralWCS(axis=coord_axis,
                               specsys='TOPOCENT',
                               ssysobs='TOPOCENT',
                               ssyssrc='TOPOCENT',
                               restfrq=None,
                               restwav=None,
                               velosys=None,
                               zsource=None,
                               velang=None,
                               bandpass_name=bandpass_name,
                               transition=None,
                               resolving_power=resolving_power)


def _update_chunk_energy_niri(chunk, header, obs_id):
    """NIRI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_niri')
    mc.check_param(chunk, Chunk)

    # No energy information is determined for darks.  The
    # latter are sometimes only identified by a 'blank' filter.  e.g.
    # NIRI 'flats' are sometimes obtained with the filter wheel blocked off.

    # The focal_plane_mask has values like f6-2pixBl (2-pixel wide blue slit
    # used with f/6 camera) or f32-6pix (6-pixel wide slit [centered]
    # with f/32 camera).   So f#-#pix with or without ‘Bl’.   But there
    # seems to be some inconsistency in the above web page and actual
    # slits.  I don’t think there are 7- or 10-pixel slits used with the
    # f32 camera.  The Gemini archive pull-down menu only gives f32-6pix
    # and f32-9pix choices.  Assume these refer to the 7 and 10 pixel
    # slits in the web page. Use the ‘disperser’ value to lookup the
    # ‘Estimated Resolving Power’ for the slit in the beam given by the
    # focal_plane_mask value.

    n_axis = 1

    filters = get_filter_name(header)
    filter_md = em.get_filter_metadata('NIRI', filters)
    filter_name = em.om.get('filter_name')

    logging.error('filters {} filter_name {}'.format(filters, filter_name))

    mode = em.om.get('mode')
    if mode == 'imaging':
        logging.debug('NIRI: SpectralWCS imaging mode for {}.'.format(obs_id))
        c_val, delta, resolving_power = _imaging_energy(filter_md)
    elif mode in ['LS', 'spectroscopy']:
        logging.debug('NIRI: SpectralWCS mode {} for {}.'.format(mode, obs_id))
        if (chunk.position is not None and chunk.position.axis is not None
            and chunk.position.axis.function is not None
            and chunk.position.axis.function.dimension is not None
                and chunk.position.axis.function.dimension.naxis1 is not None):
            n_axis = chunk.position.axis.function.dimension.naxis1
        c_val = filter_md['wl_min']  # minimum wavelength
        delta = filter_md['wl_eff_width'] / n_axis / 1.0e4
        disperser = em.om.get('disperser')
        if disperser in ['Jgrism', 'Jgrismf32', 'Hgrism', 'Hgrismf32',
                         'Kgrism', 'Kgrismf32', 'Lgrism', 'Mgrism']:
            bandpass_name = disperser[0]
            f_ratio = em.om.get('focal_plane_mask')
            logging.debug('NIRI: Bandpass name is {} f_ratio is {} for {}'.format(
                bandpass_name, f_ratio, obs_id))
            if bandpass_name in em.NIRI_RESOLVING_POWER:
                resolving_power = \
                    em.NIRI_RESOLVING_POWER[bandpass_name][f_ratio]
            else:
                logging.info('NIRI: No resolving power for {}.'.format(obs_id))
                resolving_power = None
        else:
            logging.info(
                'NIRI: Unknown disperser value {} for {}'.format(disperser,
                                                                 obs_id))
            resolving_power = None
    else:
        raise mc.CadcException(
            'NIRI: Do not understand mode {} for {}'.format(mode, obs_id))

    _build_chunk_energy(chunk, n_axis, c_val, delta, filter_name,
                        resolving_power)
    logging.debug('End _update_chunk_energy_niri')


def _update_chunk_energy_gpi(chunk, header):
    """NIRI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_gpi')
    mc.check_param(chunk, Chunk)

    n_axis = 1

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
    logging.debug('End _update_chunk_energy_gpi')


def _update_chunk_energy_f2(chunk, header):
    """NIRI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_f2')
    mc.check_param(chunk, Chunk)

    n_axis = 1
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


# select filter_id, wavelength_central, wavelength_lower, wavelength_upper
# from gsa..gsa_filters where instrument = 'NICI'
# 0 - central
# 1 - lower
# 2 - upper
#
# dict with the barcodes stripped from the names as returned by query
NICI = {'Br-gamma': [2.168600, 2.153900, 2.183300],
        'CH4-H1%L': [1.628000, 1.619300, 1.636700],
        'CH4-H1%S': [1.587000, 1.579500, 1.594500],
        'CH4-H1%Sp': [1.603000, 1.594900, 1.611100],
        'CH4-H4%L': [1.652000, 1.619000, 1.685000],
        'CH4-H4%S': [1.578000, 1.547000, 1.609000],
        'CH4-H6.5%L': [1.701000, 1.652400, 1.749600],
        'CH4-H6.5%S': [1.596000, 1.537250, 1.654750],
        'CH4-K5%L': [2.241000, 2.187500, 2.294500],
        'CH4-K5%S': [2.080000, 2.027500, 2.132500],
        'FeII': [1.644000, 1.631670, 1.656330],
        'H2-1-0-S1': [2.123900, 2.110800, 2.137000],
        'H20-Ice-L': [3.090000, 3.020000, 3.150000],
        'H': [1.650000, 1.490000, 1.780000],
        'J': [1.250000, 1.150000, 1.330000],
        'K': [2.200000, 2.030000, 2.360000],
        'Kcont': [2.271800, 2.254194, 2.289406],
        'Kprime': [2.120000, 1.950000, 2.300000],
        'Ks': [2.150000, 1.990000, 2.300000],
        'Lprime': [3.780000, 3.430000, 4.130000],
        'Mprime': [4.680000, 4.550000, 4.790000]}


def _update_chunk_energy_nici(chunk, header, obs_id):
    """NICI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_nici')
    mc.check_param(chunk, Chunk)

    n_axis = 1
    mode = em.om.get('mode')
    logging.error('mode is {}'.format(mode))

    # use filter names from headers, because there's a different
    # filter/header, and the JSON summary value obfuscates that

    temp = get_filter_name(header)
    if '+' in temp:
        raise mc.CadcException(
            'NICI: Do not understand filter {} for {}'.format(temp, obs_id))
    filter_name = temp.split('_G')[0]
    filter_md = em.get_filter_metadata('NICI', filter_name)

    if mode == 'imaging':
        logging.debug('NICI: SpectralWCS imaging mode for {}.'.format(obs_id))
        if len(filter_md) == 0:
            if filter_name in NICI:
                w_max = NICI[filter_name][2]
                w_min = NICI[filter_name][1]
                reference_wavelength = w_min
                delta = w_max - w_min
                resolving_power = (w_max + w_min)/(2*delta)
            else:
                raise mc.CadcException(
                    'NICI: Unprepared for filter {} from {}'.format(
                        filter_name, obs_id))
        else:
            reference_wavelength, delta, resolving_power = \
                _imaging_energy(filter_md)
    else:
        raise mc.CadcException(
            'NICI: Do not understand mode {} from {}'.format(mode, obs_id))

    # NICI has two different bandpass names (most of the time) in its two
    # chunks.  Pat says in this case nothing will be put in the bandpass
    # name for the plane.  Add code to combine the two chunk bandpass names
    # to create a plane bandpass name only for this instrument
    #
    # DB - 14-02-19 value looks clearer (as two filters) with a space on
    # both sides of the ‘+’.
    
    temp = em.om.get('filter_name')
    _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                        temp.replace('+', ' + '), resolving_power)

    logging.debug('End _update_chunk_energy_nici')


def _update_chunk_energy_trecs(chunk, header, obs_id):
    """TReCS-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_trecs')
    mc.check_param(chunk, Chunk)

    # Look at the json ‘disperser’ value.  If = LowRes-20" then
    # resolving power = 80.  If LowRes-10" then resolving poewr = 100.
    # There was a high-res mode but perhaps never used.  Again,
    # naxis = NAXIS1 header value and other wcs info as above for NIRI.  But
    # might take some string manipulation to match filter names with SVO
    # filters.

    orig_filter_name = em.om.get('filter_name')
    filter_name = orig_filter_name
    temp = filter_name.split('-')
    if len(temp) > 1:
        filter_name = temp[0]
    # what filter names look like in Gemini metadata, and what they
    # look like at the SVO, aren't necessarily the same
    trecs_repair = {'K': 'k',
                    'L': 'l',
                    'M': 'm',
                    'N': 'n',
                    'Nprime': 'nprime',
                    'Qw': 'Qwide',
                    'NeII_ref2': 'NeII_ref'}
    if filter_name in trecs_repair:
        filter_name = trecs_repair[filter_name]
    logging.debug('TReCS: filter_name is {} for {}'.format(filter_name, obs_id))
    filter_md = em.get_filter_metadata('TReCS', filter_name)

    resolving_power = None
    mode = em.om.get('mode')
    if mode in ['imaging', 'IFP', 'IFS']:
        logging.debug('TRECS: imaging mode for {}.'.format(obs_id))
        n_axis = 1
        reference_wavelength, delta, resolving_power = \
            _imaging_energy(filter_md)
    elif mode in ['LS', 'spectroscopy', 'MOS']:
        logging.debug('TReCS: LS|Spectroscopy mode for {}.'.format(obs_id))
        n_axis = header.get('NAXIS1')
        delta = filter_md['wl_eff_width'] / n_axis / 1.0e4
        reference_wavelength = em.om.get('central_wavelength')
        disperser = em.om.get('disperser')
        if disperser is not None:
            if disperser == 'LowRes-20':
                resolving_power = 80.0
            elif disperser == 'LowRes-10':
                resolving_power = 100.0
    else:
        raise mc.CadcException(
            'Do not understand mode {} for {}'.format(mode, obs_id))

    _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                        orig_filter_name, resolving_power)
    logging.debug('End _update_chunk_energy_trecs')


def _update_chunk_energy_michelle(chunk, obs_id):
    """Michelle-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_michelle')
    mc.check_param(chunk, Chunk)

    # what filter names look like in Gemini metadata, and what they
    # look like at the SVO, aren't necessarily the same

    n_axis = 1
    orig_filter_name = em.om.get('filter_name')
    filter_name = orig_filter_name
    temp = filter_name.split('-')
    if len(temp) > 1:
        filter_name = temp[0]

    michelle_repair = {'I79B10': 'Si1',
                       'I88B10': 'Si2',
                       'I97B10': 'Si3',
                       'I103B10': 'Si4',
                       'I105B53': 'N',
                       'I112B21': 'Np',
                       'I116B9': 'Si5',
                       'I125B9': 'Si6',
                       'I185B9': 'Qa',
                       'I209B42': 'Q'}
    if filter_name in michelle_repair:
        filter_name = michelle_repair[filter_name]
    filter_md = em.get_filter_metadata('michelle', filter_name)

    mode = em.om.get('mode')
    logging.debug('michelle: mode is {}, filter_name is {}, for {}'.format(
        mode, filter_name, obs_id))
    if mode in ['imaging', 'IFP', 'IFS']:
        logging.debug(
            'michelle: SpectralWCS imaging mode for {}.'.format(obs_id))
        reference_wavelength, delta, resolving_power = \
            _imaging_energy(filter_md)
    elif mode in ['LS', 'spectroscopy', 'MOS']:
        logging.debug(
            'michelle: SpectralWCS LS|Spectroscopy mode for {}.'.format(
                obs_id))
        delta = filter_md['wl_eff_width'] / n_axis / 1.0e4
        reference_wavelength = em.om.get('central_wavelength')
    else:
        raise mc.CadcException(
            'michelle: Do not understand mode {} for {}'.format(mode, obs_id))

    _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                        orig_filter_name, resolving_power=None)
    logging.debug('End _update_chunk_energy_michelle')


def _update_chunk_energy_nifs(chunk, header, obs_id):
    """NIFS-specific chunk-level Spectral WCS construction."""
    logging.debug('Begin _update_chunk_energy_nifs')
    mc.check_param(chunk, Chunk)

    # Imaging: treated as for other instruments.
    #
    # Spectroscopy:
    #
    # Need grating name from json ‘disperser’ value
    # Need central wavelength from json ‘central_wavelength’ value
    # Need filter name from json ‘filter’ value
    # Need header NAXIS1 value in extension for number of pixels in
    # dispersion direction
    #
    # Then the top table on this page to create a lookup for
    # upper/lower wavelength values and spectral resolution for
    # the given grating/filter combination:
    # https://www.gemini.edu/sciops/instruments/nifs/ifu-spectroscopy/gratings.
    #
    # cdelta will be as for other instruments (upper - lower)/naxis

    # Use the filter names WITHOUT the ‘_298K’ suffix (since those are for
    # the cold filters at operating temperature).

    # key = grating name
    # key = associated filter
    # 0 = Central Wavelength (microns),
    # 1, 2 = Spectral Range,
    # 3 = Spectral Resolution,
    # 4 = Velocity Resolution (km/s)

    # From the second table:
    # Table 2: The NIFS gratings can be tuned to different central
    # wavelengths. The short and long limits of the possible tuned
    # central wavelengths and the associated filters required are:
    # Grating, Name, Short Wavelength Limit (microns), Long Wavelength Limit (microns)	Short Wavelength Filter	Long Wavelength Filter
    # Z	0.94	1.16	ZJ	ZJ
    # J	1.14	1.36	ZJ	JH
    # H	1.48	1.82	JH	HK
    # K	1.98	2.41	HK	HK

    nifs_lookup = {'Z':       {'ZJ': [1.05, 0.94, 1.15, 4990.0, 60.1]},
                   'J':       {'ZJ': [1.25, 1.15, 1.33, 6040.0, 49.6]},
                   'H':       {'JH': [1.65, 1.49, 1.80, 5290.0, 56.8]},
                   'K':       {'HK': [2.20, 1.99, 2.40, 5290.0, 56.7]},
                   'K_Short': {'HK': [2.20, 1.98, 2.41, None, None]},
                   'K_Long':  {'HK': [2.20, 1.98, 2.41, None, None]}}

    filter_name = em.om.get('filter_name')
    filter_md = em.get_filter_metadata('NIFS', filter_name)

    spectroscopy = em.om.get('spectroscopy')
    logging.debug('NIFS: spectroscopy is {}, filter is {} for {}'.format(
        spectroscopy, filter_name, obs_id))
    if spectroscopy is not None:
        if spectroscopy:
            logging.debug('NIFS: spectroscopy for {}.'.format(obs_id))
            n_axis = header.get('NAXIS1')
            reference_wavelength = em.om.get('central_wavelength')
            grating = em.om.get('disperser')
            if grating in nifs_lookup:
                if filter_name in nifs_lookup[grating]:
                    wl_min = nifs_lookup[grating][filter_name][1]
                    wl_max = nifs_lookup[grating][filter_name][2]
                    delta = (wl_max - wl_min) / n_axis
                    resolving_power = nifs_lookup[grating][filter_name][3]
                else:
                    raise mc.CadcException('NIFS: mystery filter_name {} for {}'.format(
                        filter_name, obs_id))
            else:
                raise mc.CadcException(
                    'NIFS: mystery grating {} for {}'.format(grating, obs_id))
        else:
            logging.debug('NIFS: imaging for {}.'.format(obs_id))
            n_axis = 1
            reference_wavelength, delta, resolving_power = \
                _imaging_energy(filter_md)
    else:
        raise mc.CadcException(
            'NIFS: No spectroscopy information for {}'.format(obs_id))

    _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                        filter_name, resolving_power)
    logging.debug('End _update_chunk_energy_nifs')


def _update_chunk_energy_gnirs(chunk, header):
    """GNIRS-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_gnirs')
    mc.check_param(chunk, Chunk)

    n_axis = header.get('NAXIS2')
    filter_name = em.om.get('filter_name')
    spectroscopy = em.om.get('spectroscopy')
    # spectroscopy mode
    if spectroscopy:
        logging.debug('gnirs: SpectralWCS Spectroscopy mode.')
        reference_wavelength = em.om.get('central_wavelength')
        disperser = em.om.get('disperser')
        if 'XD' in disperser:
            logging.debug('gnirs: cross dispersed mode.')
            # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/
            # crossdispersed-xd-spectroscopy/xd-prisms
            xd_mode = {'SB+SXD': [0.9, 2.5],
                       'LB+LXD': [0.9, 2.5],
                       'LB+SXD': [1.2, 2.5]}
            bounds = None
            coverage = disperser[-3:]
            camera = em.om.get('camera')
            if camera.startswith('Short'):
                bounds = xd_mode['{}+{}'.format('SB', coverage)]
            elif camera.startswith('Long'):
                bounds = xd_mode['{}+{}'.format('LB', coverage)]
            if not bounds:
                raise mc.CadcException(
                    'gnirs: Do not understand xd mode {} {}'
                        .format(camera, coverage))

        else:
            logging.debug('gnirs: long slit mode.')
            # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy
            long_slit_mode = {'11': {'X': [1.03, 1.17],
                                     'J': [1.17, 1.37],
                                     'H': [1.47, 1.80],
                                     'K': [1.91, 2.49],
                                     'L': [2.80, 4.20],
                                     'M': [4.40, 6.00]},
                              '32': {'X': [1.03, 1.17],
                                     'J': [1.17, 1.37],
                                     'H': [1.49, 1.80],
                                     'K': [1.91, 2.49],
                                     'L': [2.80, 4.20],
                                     'M': [4.40, 6.00]},
                              '111': {'X': [1.03, 1.17],
                                      'J': [1.17, 1.37],
                                      'H': [1.49, 1.80],
                                      'K': [1.91, 2.49],
                                      'L': [2.80, 4.20],
                                      'M': [4.40, 6.00]}}
            bandpass = filter_name[0]
            grating = disperser.split('_')[0]
            bounds = long_slit_mode[grating][bandpass]
            if not bounds:
                raise mc.CadcException(
                    'gnirs: Do not understand long slit mode with {} {}'
                        .format(bandpass, bounds))

        delta = (bounds[1] - bounds[0]) / n_axis

    # imaging mode
    else:
        logging.debug('gnirs: SpectralWCS imaging mode.')
        # https://www.gemini.edu/sciops/instruments/gnirs/imaging
        imaging = {'Y': [0.97, 1.07],
                   'J': [1.17, 1.33],
                   'J order blocking)': [1.17, 1.37],
                   'H': [1.49, 1.80],
                   'K': [2.03, 2.37],
                   'K order blocking': [1.90, 2.49],
                   'H2': [2.105, 2.137],
                   'PAH': [3.27, 3.32]}

        if len(filter_name) == 1 or filter_name is 'H2' or filter_name is 'PAH':
            bandpass = filter_name
        else:
            bandpass = filter_name[0]

        bounds = imaging[bandpass]

        reference_wavelength = bounds[0]
        delta = (bounds[1] - bounds[0]) / n_axis

    _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                        filter_name, resolving_power=None)
    logging.debug('End _update_chunk_energy_gnirs')


# select filter_id, wavelength_central, wavelength_lower, wavelength_upper
# from gsa..gsa_filters where instrument = 'NICI'
# 0 - central
# 1 - lower
# 2 - upper
#
# dict with the filter wheels stripped from the names as returned by query
PHOENIX = {'2030': [4.929000, 4.808000, 5.050000],
           '2150': [4.658500, 4.566000, 4.751000],
           '2462': [4.078500, 4.008000, 4.149000],
           '2734': [3.670500, 3.610000, 3.731000],
           '2870': [3.490500, 3.436000, 3.545000],
           '3010': [3.334500, 3.279000, 3.390000],
           '3100': [3.240000, 3.180000, 3.300000],
           '3290': [3.032500, 2.980000, 3.085000],
           '4220': [2.370000, 2.348000, 2.392000],
           '4308': [2.322500, 2.296000, 2.349000],
           '4396': [2.272500, 2.249000, 2.296000],
           '4484': [2.230000, 2.210000, 2.250000],
           '4578': [2.185000, 2.160000, 2.210000],
           '4667': [2.143000, 2.120000, 2.166000],
           '4748': [2.104000, 2.082000, 2.126000],
           '6073': [1.647000, 1.632000, 1.662000],
           '6420': [1.557500, 1.547000, 1.568000],
           '7799': [1.280500, 1.271000, 1.290000],
           '8265': [1.204500, 1.196000, 1.213000],
           '9232': [1.083000, 1.077000, 1.089000],
           'L2870': [3.490500, 3.436000, 3.545000]}


def _update_chunk_energy_phoenix(chunk, header, obs_id):
    """Phoenix-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_phoenix')
    mc.check_param(chunk, Chunk)

    # I forgot what we used for NAXIS for TReCS (i.e. NAXIS1 or 2
    # value).  Looks like you’ll have to look at the header FILT_POS
    # value to determine the filter since it’s not in json metadata.
    # DB - 12-02-19 - Phoenix should be the same but uses NAXIS2 for the
    # length of the dispersion axis.

    # DB - 12-02-19 - Note that the parenthetical numbers
    # after the Phoenix filter names (in the header) indicates the
    # filter wheel slot the filter is in and may occasionally change
    # so should be disregarded.
    filter_name = get_filter_name(header, 'FILT_POS')
    if len(filter_name) > 0:
        filter_name = filter_name.split()[0]

    logging.debug(
        'Phoenix: filter_name is {} for {}'.format(filter_name, obs_id))

    resolving_power = None
    spectroscopy = em.om.get('spectroscopy')
    if spectroscopy is None:
        raise mc.CadcException(
            'Phoenix: Spectroscopy undefined for {}'.format(obs_id))
    else:
        if spectroscopy:
            logging.debug(
                'Phoenix: LS|Spectroscopy mode for {}.'.format(obs_id))
            global phoenix_spectral_axis
            n_axis = phoenix_spectral_axis
            if filter_name in PHOENIX:
                w_max = PHOENIX[filter_name][2]
                w_min = PHOENIX[filter_name][1]
                delta = (w_min + w_max) / n_axis / 1.0e4
                reference_wavelength = PHOENIX[filter_name][0]
            elif len(filter_name) == 0:
                w_max = 100000.0
                w_min = 0.0
                delta = (w_min + w_max) / n_axis / 1.0e4
                reference_wavelength = (w_max + w_min)/2.0
            else:
                raise mc.CadcException(
                    'Phoenix: mystery filter name {} for {}'.format(
                        filter_name, obs_id))
        else:
            logging.debug('Phoenix: imaging mode for {}.'.format(obs_id))
            n_axis = 1
            if filter_name in PHOENIX:
                w_max = PHOENIX[filter_name][2]
                w_min = PHOENIX[filter_name][1]
            elif len(filter_name) == 0:
                w_max = 100000.0
                w_min = 0.0
            else:
                raise mc.CadcException(
                    'Phoenix: mystery filter name {} for {}'.format(
                        filter_name, obs_id))
            filter_md = {'wl_min': w_min,
                         'wl_eff_width': w_max - w_min}
            reference_wavelength, delta, resolving_power = \
                _imaging_energy(filter_md)

    _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                        filter_name, resolving_power)
    logging.debug('End _update_chunk_energy_phoenix')


def _update_chunk_energy_flamingos(chunk, header, data_product_type, obs_id):
    """Flamingos-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_flamingos')
    mc.check_param(chunk, Chunk)

    # DB - 18-02-19 - FLAMINGOS spectral WCS should be similar to what was
    # done for NIRI.  Use NAXIS1 keyword value to determine number of
    # pixels.  The GRISM and FILTER keywords give the same filter ID.
    # The SVO filter information for KPNO/Flamingos has a blue leak in the
    # JK blocking filter which give too large a spectral range.  Can you
    # instead hard code min/max wavelengths using the top two lines in
    # this table 7:
    # http://www-kpno.kpno.noao.edu/manuals/flmn/flmn.user.html#flamspec.
    # Use these for min/max wavelengths and use the average as the ‘central’
    # wavelength.  Spectral resolution is fixed at about 1300 for both grisms.

    # spectral wcs units are microns, values from Table 7 are angstroms
    lookup = {'JH': [8910, 18514],
              'HK': [10347, 27588]}

    filter_name = get_filter_name(header, 'FILTER')
    if filter_name in ['JH', 'HK']:
        filter_md = {'wl_min': lookup[filter_name][0] / 1.0e4,
                     'wl_eff_width': (lookup[filter_name][1] -
                                      lookup[filter_name][0]) / 1.0e4}
        resolving_power = 1300.0 / 1.0e4
    else:
        filter_md = em.get_filter_metadata('Flamingos', filter_name)
        resolving_power = None

    if data_product_type == DataProductType.SPECTRUM:
        logging.debug('Flamingos: SpectralWCS for {}.'.format(obs_id))
        n_axis = header.get('NAXIS1')
        reference_wavelength = float(filter_md['wl_min'])  # minimum wavelength
        delta = filter_md['wl_eff_width'] / n_axis
    elif data_product_type == DataProductType.IMAGE:
        logging.debug(
            'Flamingos: SpectralWCS imaging mode for {}.'.format(obs_id))
        n_axis = 1
        reference_wavelength, delta, resolving_power = \
            _imaging_energy(filter_md)
    else:
        raise mc.CadcException(
            'Flamingos: mystery data product type {} for {}'.format(
                data_product_type, obs_id))

    _build_chunk_energy(chunk, n_axis, reference_wavelength, delta,
                        filter_name, resolving_power)
    logging.debug('End _update_chunk_energy_flamingos')


def _reset_energy(data_label):
    """
    Return True if there should be no energy WCS information created at
    the chunk level.

    :param data_label str for useful logging information only.
    """
    result = False
    observation_type = em.om.get('observation_type')
    om_filter_name = em.om.get('filter_name')

    if ((observation_type is not None and observation_type == 'DARK') or
        (om_filter_name is not None and ('blank' in om_filter_name or
         'Blank' in om_filter_name))):
        logging.info(
            'No chunk energy for {} obs type {} filter name {}'.format(
                data_label, observation_type, om_filter_name))
        result = True
    return result


def _imaging_energy(filter_md):
    """"""
    # all the filter units are Angstroms, and the chunk-level energy WCS
    # unit of choice is 'um', therefore 1.0e4

    # NAXIS = 1
    # CTYPE = ‘WAVE’
    # CUNIT = ‘um’
    # CRPIX = NAXIS / 2.0
    # CRVAL = wl_min
    # CDELT = wl_max - wl_min
    # resolving power = (wl_max  + wl_min)/(2*CDELT)
    # bandpass_name = filter name
    # CDELT is the SVO’s effective width, W_effective in their tables which
    # corresponds more or less to wl_max - wl_min.

    # CRPIX values should all be 0.5 for imaging spectral WCS as long as
    # you use the lower wavelength boundary of the filter for the
    # corresponding CRVAL. DB - 02-13-19

    c_val = filter_md['wl_min'] / 1.0e4
    delta = filter_md['wl_eff_width'] / 1.0e4
    resolving_power = c_val / delta
    return c_val, delta, resolving_power


def get_filter_name(header, lookup='FILTER'):
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
        if lookup in key:
            value = header.get(key).lower()
            ignore = False
            for ii in filters2ignore:
                if ii.startswith(value) or value.startswith(ii):
                    ignore = True
                    break
            if ignore:
                continue
            else:
                header_filters.append(header.get(key).strip())
        filters = '+'.join(header_filters)
    logging.info('Filters are {}'.format(filters))
    return filters


def _update_position_from_zeroth_header(artifact, headers, log_id):
    """Make the 0th header spatial WCS the WCS for all the
    chunks."""
    primary_header = headers[0]

    # naxis values are only available from extensions

    primary_header['NAXIS1'] = headers[-1].get('NAXIS1')
    primary_header['NAXIS2'] = headers[-1].get('NAXIS2')

    wcs_parser = WcsParser(primary_header, 'file name', 0)
    primary_chunk = Chunk()
    wcs_parser.augment_position(primary_chunk)

    for part in artifact.parts:
        if part == '0':
            continue
        for chunk in artifact.parts[part].chunks:
            if primary_chunk.position is not None:
                chunk.position = primary_chunk.position
                chunk.position_axis_1 = 1
                chunk.position_axis_2 = 2


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
