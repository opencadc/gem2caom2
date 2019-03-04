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
import math
import os
import sys
import re
import traceback

from astropy import units
from astropy.coordinates import SkyCoord

from caom2 import Observation, ObservationIntentType, DataProductType
from caom2 import CalibrationLevel, TargetType, ProductType, Chunk, Axis
from caom2 import SpectralWCS, CoordAxis1D, CoordFunction1D, RefCoord
from caom2 import TypedList, CoordRange1D
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2utils import WcsParser
from caom2pipe import manage_composable as mc
from caom2pipe import astro_composable as ac

import gem2caom2.external_metadata as em
from gem2caom2.gem_name import GemName, COLLECTION, ARCHIVE, SCHEME
from gem2caom2.svofps import FilterMetadata

__all__ = ['main_app2', 'update', 'APPLICATION']

# TODO LIST:
#
# PI information
# spectroscopy check

APPLICATION = 'gem2caom2'


# GPI radius == 2.8 arcseconds, according to Christian Marois via DB 02-07-19
# DB - 18-02-19 - Replace “5.0” for “2.8" for GPI field of view.
# NOTE:  To be more accurate for GRACES this size could be reduced
# from 5" to 1.2" since that’s the size of the fibre.
# OSCIR - http://www.gemini.edu/sciops/instruments/oscir/oscirIndex.html
# bHROS - DB - 20-02-19 - bHROS ‘bounding box’ is only 0.9".
#                         A very small fibre.
# HOKUPAA - http://www.gemini.edu/sciops/instruments/uhaos/uhaosIndex.html
RADIUS_LOOKUP = {em.Inst.GPI: 2.8 / 3600.0,  # units are arcseconds
                 em.Inst.GRACES: 1.2 / 3600.0,
                 em.Inst.PHOENIX: 5.0 / 3600.0,
                 em.Inst.OSCIR: 11.0 / 3600.0,
                 em.Inst.HOKUPAA: 4.0 / 3600.0,
                 em.Inst.BHROS: 0.9 / 3600.0}


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
            instrument = em.Inst(header.get('INSTRUME'))
            if instrument is not None:
                if instrument == em.Inst.PHOENIX:
                    if _is_phoenix_calibration(header):
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                elif instrument == em.Inst.OSCIR:
                    if _is_oscir_calibration(header):
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                elif instrument == em.Inst.FLAMINGOS:
                    if _is_flamingos_calibration():
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                elif instrument == em.Inst.HOKUPAA:
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
    instrument = _get_instrument()
    if instrument == em.Inst.HOKUPAA:
        result = _get_pix_scale(header)
    elif instrument == em.Inst.OSCIR:
        result = 0.0890 / 3600.0
    elif instrument == em.Inst.GPI:
        # FOV is 2.8" x 2.8" on each side or 2.8/3600.0 degrees
        # cd1_1 = 2.8/(3600.0 * naxis1)
        result = 2.8 / (3600.0 * header.get('NAXIS1'))
    else:
        cdelt1 = header.get('CDELT1')
        if cdelt1 is None:
            result = RADIUS_LOOKUP[instrument]
        else:
            result = cdelt1
    return result


def get_cd22(header):
    instrument = _get_instrument()
    if instrument == em.Inst.HOKUPAA:
        result = _get_pix_scale(header)
    elif instrument == em.Inst.OSCIR:
        result = 0.0890 / 3600.0
    elif instrument == em.Inst.GPI:
        # FOV is 2.8" x 2.8" on each side or 2.8/3600.0 degrees
        # cd2_2 = 2.8/(3600.0 * naxis2)
        result = 2.8 / (3600.0 * header.get('NAXIS2'))
    else:
        cdelt2 = header.get('CDELT2')
        if cdelt2 is None:
            result = RADIUS_LOOKUP[instrument]
        else:
            result = cdelt2
    return result


def get_crpix1(header):
    instrument = _get_instrument()
    if instrument == em.Inst.HOKUPAA:
        crpix1 = header.get('CRPIX1')
        if crpix1 is None:
            result = None
        else:
            result = crpix1 / 2.0
    elif instrument in [em.Inst.OSCIR, em.Inst.GPI]:
        naxis1 = header.get('NAXIS1')
        if naxis1 is None:
            result = None
        else:
            result = naxis1 / 2.0
    else:
        result = 1.0
    return result


def get_crpix2(header):
    instrument = _get_instrument()
    if instrument == em.Inst.HOKUPAA:
        crpix2 = header.get('CRPIX2')
        if crpix2 is None:
            result = None
        else:
            result = crpix2 / 2.0
    elif instrument in [em.Inst.OSCIR, em.Inst.GPI]:
        naxis2 = header.get('NAXIS2')
        if naxis2 is None:
            result = None
        else:
            result = naxis2 / 2.0
    else:
        result = 1.0
    return result


def get_data_product_type(header):
    """
    Calculate the Plane DataProductType.

    :param header:  The FITS header for the current extension.
    :return: The Plane DataProductType, or None if not found.
    """
    instrument = _get_instrument()
    if instrument == em.Inst.FLAMINGOS:
        result, ignore = _get_flamingos_mode(header)
    else:
        mode = em.om.get('mode')
        obs_type = _get_obs_type(header)
        if mode is None:
            raise mc.CadcException('No mode information found for {}'.format(
                em.om.get('filename')))
        elif instrument == em.Inst.GPI:
            logging.error('mode is {}'.format(mode))
            # DB - 22-02-19 FOR GPI only:  To determine if the data type
            # is an ‘image’ or ‘spectrum’:
            #     json ‘mode’ = IFP then ‘image’
            #     json ‘mode’ = IFS then ‘spectrum’
            if mode in ['IFP', 'imaging']:
                result = DataProductType.IMAGE
            elif mode == 'IFS':
                result = DataProductType.SPECTRUM
            else:
                raise mc.CadcException(
                    'Mystery GPI mode {} for {}'.format(
                        mode, em.om.get('filename')))
        elif ((mode == 'imaging') or
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


def get_dec(header):
    """
    Get the declination. Rely on the JSON metadata, because it's all in
    the same units (degrees).

    :param header:  The FITS header for the current extension.
    :return: declination, or None if not found.
    """
    instrument = _get_instrument()
    if instrument == em.Inst.HOKUPAA:
        ra, dec = _get_sky_coord(header, 'RA', 'DEC')
        result = dec
    elif instrument == em.Inst.OSCIR:
        ra, dec = _get_sky_coord(header, 'RA_TEL', 'DEC_TEL')
        result = dec
    elif instrument == em.Inst.BHROS:  # ra/dec not in json
        result = header.get('DEC')
    else:
        result = em.om.get('dec')
    return result


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
    result = em.om.get('exposure_time')
    instrument = _get_instrument()
    if instrument == em.Inst.OSCIR:
        # DB - 20-02-19 - json ‘exposure_time’ is in minutes, so multiply
        # by 60.0.
        result = result * 60.0
    elif instrument == em.Inst.FLAMINGOS:
        # exposure_time is null in the JSON
        result = mc.to_float(header.get('EXP_TIME'))
    return result


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
    cal_values = ['GCALflat', 'Bias', 'BIAS', 'Twilight', 'Ar', 'FLAT',
                  'flat', 'ARC', 'Domeflat', 'DARK', 'dark']
    dl = header.get('DATALAB')
    lookup = _get_obs_class(header)
    logging.debug('observation_class is {} for {}'.format(lookup, dl))
    if lookup is None:
        type_lookup = _get_obs_type(header)
        logging.debug('observation_type is {} for {}'.format(type_lookup, dl))
        if type_lookup is None:
            data_label = em.om.get('data_label')
            logging.debug('data_label is {}'.format(data_label))
            if (data_label is None or
                    (data_label is not None and '-CAL' not in data_label)):
                object_value = header.get('OBJECT')
                logging.debug(
                    'object_value is {} for {}'.format(object_value, dl))
                if object_value is not None:
                    if object_value in cal_values:
                        result = ObservationIntentType.CALIBRATION
                    else:
                        # check that an individual cal value is not part of
                        # the object_value
                        cal_value_found = False
                        for ii in cal_values:
                            if ii in object_value:
                                result = ObservationIntentType.CALIBRATION
                                cal_value_found = True
                                break
                        if not cal_value_found:
                            result = ObservationIntentType.SCIENCE
            else:
                if '-CAL' in data_label:
                    result = ObservationIntentType.CALIBRATION
                else:
                    result = ObservationIntentType.SCIENCE
        else:
            if type_lookup in cal_values:
                result = ObservationIntentType.CALIBRATION
            else:
                instrument = em.Inst(header.get('INSTRUME'))
                logging.error('instrument is {}'.format(instrument))
                if instrument is not None:
                    if instrument == em.Inst.TRECS:
                        data_label = header.get('DATALAB')
                        if data_label is not None and '-CAL' in data_label:
                            result = ObservationIntentType.CALIBRATION
                    else:
                        result = ObservationIntentType.SCIENCE
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
    instrument = _get_instrument()
    if instrument == em.Inst.PHOENIX:
        result = _get_phoenix_obs_type(header)
    elif instrument == em.Inst.HOKUPAA:
        # DB - 01-18-19 Use IMAGETYP as the keyword
        result = header.get('IMAGETYP')
    return result


def get_ra(header):
    """
    Get the right ascension. Rely on the JSON metadata, because it's all in
    the same units (degrees).

    :param header:  The FITS header for the current extension.
    :return: ra, or None if not found.
    """
    instrument = _get_instrument()
    if instrument == em.Inst.HOKUPAA:
        ra, dec = _get_sky_coord(header, 'RA', 'DEC')
        result = ra
    elif instrument == em.Inst.OSCIR:
        ra, dec = _get_sky_coord(header, 'RA_TEL', 'DEC_TEL')
        result = ra
    elif instrument == em.Inst.BHROS:
        result = header.get('RA')  # ra/dec not in json
        logging.error('getting ra from here {}'.format(result))
    else:
        result = em.om.get('ra')
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
    instrument = _get_instrument()
    if instrument in [em.Inst.FLAMINGOS, em.Inst.OSCIR]:
        # Another FLAMINGOS correction needed:  DATE-OBS in header and
        # json doesn’t include the time  but sets it to 00:00:00.  You have
        # two choices:  concatenate DATE-OBS and UTC header values to the
        # standard form “2002-11-06T07:06:00.3” and use as you use json
        # value for computing temporal WCS data or, for FLAMINGOS only,
        # use the MJD header value in the header as the CRVAL for time.

        # DB - 20-02-19 - OSCIR json ‘ut_datetime’ is not correct.  Must
        # concatenate DATE-OBS and UTC1 values and convert to MJD as usual
        # or use MJD directly (seems to be correct starting MJD)

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


def _get_pix_scale(header):
    pix_scale = header.get('PIXSCALE')
    if pix_scale is None:
        result = None
    else:
        result = pix_scale / 3600.0
    return result


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
        else:
            # DB - 27-02-19 - Default to 'OBJECT'
            obs_type = 'OBJECT'
    return data_type, obs_type


def _get_instrument():
    return em.Inst(em.om.get('instrument'))


def _get_sky_coord(header, ra_key, dec_key):
    ra_hours = header.get(ra_key)
    dec_hours = header.get(dec_key)
    if ra_hours is None or dec_hours is None:
        ra_deg = None
        dec_deg = None
    else:
        result = SkyCoord('{} {}'.format(ra_hours, dec_hours),
                          unit=(units.hourangle, units.deg))
        ra_deg = result.ra.degree
        dec_deg = result.dec.degree
    return ra_deg, dec_deg


def _get_phoenix_obs_type(header):
    # DB - 18-02-19 - make PHOENIX searches more useful by making estimated
    # guesses at observation type
    #
    # Looking at a few random nights it looks like a reasonable way to
    # determine if an observation is a FLAT is to look for ‘flat’ (any
    # case combination) in the json ‘object’ value or ‘gcal’ or ‘GCAL’.
    # But if ‘gcal’ AND ‘dark’ are in the ‘object’ string it’s a DARK.
    # (see below)

    # Ditto for an ARC if json ‘object’ contains the string ‘arc’

    # It’s a DARK obs type if the value of FITS header VIEW_POS contains
    # the string ‘dark’ OR if ‘dark’ is in json ‘object’ string.
    # (There appear to be cases where darks are taken with the
    # VIEW_POS = open.)  When it’s a dark do NOT generate WCS for energy
    # since the filter is often ‘open’ and energy info isn’t important for
    # DARK exposures.

    result = 'OBJECT'
    object_value = em.om.get('object').lower()
    view_pos = header.get('VIEW_POS')
    if 'flat' in object_value:
        result = 'FLAT'
    elif ('dark' in object_value or
          (view_pos is not None and 'dark' in view_pos)):
        result = 'DARK'
    elif 'gcal' in object_value:
        result = 'FLAT'
    elif 'arc' in object_value:
        result = 'ARC'
    return result


def _is_flamingos_calibration():
    # DB - 28-02-19
    # Another relatively minor thing for FLAMINGOS:  get_obs_intent likely
    # has to have FLAMINGOS added to it.  If obs_type is DARK, ACQUISIITON,
    # FLAT or ARC then the intent is calibration, otherwise science.
    # cal_values should have those obs_type values in it as well.
    # So if ‘CAL’ is not in the data_label and obs_type is not in cal_values
    # then it’s science.
    object_value = em.om.get('object')
    for ii in ['Dark', 'DARK', 'Arc', 'flat', 'Flat', 'Acq', 'acq']:
        if ii in object_value:
            return True
    data_label = em.om.get('data_label')
    if data_label is not None and 'CAL' in data_label:
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

    instrument = _get_instrument()
    mode = em.om.get('mode')
    if not (instrument in [em.Inst.GPI, em.Inst.PHOENIX, em.Inst.HOKUPAA,
                           em.Inst.OSCIR, em.Inst.BHROS] or
            (instrument == em.Inst.GRACES and mode is not None and
                mode != 'imaging')):
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
    bp.set('Chunk.time.axis.function.refCoord.val',
           'get_time_function_val(header)')

    logging.debug('Done accumulate_fits_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes, must
    have this signature for import_module loading and execution.

    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.debug('Begin update.')
    mc.check_param(observation, Observation)

    headers = None
    if 'headers' in kwargs:
        headers = kwargs['headers']

    try:
        for p in observation.planes:
            for a in observation.planes[p].artifacts:

                instrument = em.Inst(observation.instrument.name)
                if (instrument in
                        [em.Inst.MICHELLE, em.Inst.TRECS, em.Inst.NIFS,
                         em.Inst.GNIRS]):
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

                    if part == '2' and instrument == em.Inst.GPI:
                        # GPI data sets have two extensions. First is science
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
                        if _reset_energy(observation.type, p, instrument):
                            c.energy = None
                            c.energy_axis = None
                        else:
                            filter_name = get_filter_name(
                                headers[0], header, observation.observation_id,
                                instrument)
                            if instrument == em.Inst.NIRI:
                                _update_chunk_energy_niri(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument == em.Inst.GPI:
                                _update_chunk_energy_gpi(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument == em.Inst.F2:
                                _update_chunk_energy_f2(
                                    c, headers[0],
                                    observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument in [em.Inst.GSAOI, em.Inst.HRWFS]:
                                _update_chunk_energy_general(
                                    c, instrument, [DataProductType.IMAGE],
                                    observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument == em.Inst.NICI:
                                _update_chunk_energy_nici(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument == em.Inst.TRECS:
                                _update_chunk_energy_trecs(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument == em.Inst.MICHELLE:
                                _update_chunk_energy_general(
                                    c, instrument,
                                    [DataProductType.IMAGE,
                                     DataProductType.SPECTRUM],
                                    observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument == em.Inst.NIFS:
                                _update_chunk_energy_nifs(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument == em.Inst.GNIRS:
                                _update_chunk_energy_gnirs(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument == em.Inst.PHOENIX:
                                _update_chunk_energy_phoenix(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument == em.Inst.FLAMINGOS:
                                _update_chunk_energy_flamingos(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                                ignore, obs_type = _get_flamingos_mode(header)
                                if (obs_type is not None and
                                        observation.type is None):
                                    observation.type = obs_type
                            elif instrument == em.Inst.HOKUPAA:
                                _update_chunk_energy_hokupaa(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument == em.Inst.OSCIR:
                                _update_chunk_energy_oscir(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument == em.Inst.BHROS:
                                _update_chunk_energy_bhros(
                                    c, header,
                                    observation.planes[p].data_product_type,
                                    observation.observation_id)
                            elif instrument == em.Inst.GRACES:
                                _update_chunk_energy_graces(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id)
                            elif instrument in [em.Inst.GMOS, em.Inst.GMOSN,
                                                em.Inst.GMOSS]:
                                _update_chunk_energy_gmos(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name,
                                    instrument)

                        # position WCS
                        mode = em.om.get('mode')
                        if (instrument in [em.Inst.PHOENIX, em.Inst.HOKUPAA,
                                           em.Inst.OSCIR] or
                                (instrument == em.Inst.GRACES and
                                 mode is not None and mode != 'imaging')):
                            _update_chunk_position(
                                c, header, instrument, int(part),
                                observation.observation_id)
                        elif instrument == em.Inst.BHROS:
                            _update_chunk_position(
                                c, headers[0], instrument,
                                int(part), observation.observation_id)
                        elif instrument == em.Inst.GPI:
                            _update_chunk_position(
                                c, headers[1], instrument, int(part),
                                observation.observation_id)
                            if part == '1':
                                # equinox information only available from
                                # 0th header
                                c.position.equinox = headers[0].get('TRKEQUIN')
                        elif instrument == em.Inst.FLAMINGOS:
                            _update_chunk_position_flamingos(
                                c, header, observation.observation_id)

                        # time WCS
                        if instrument == em.Inst.F2:
                            _update_chunk_time_f2(c, observation.observation_id)

        program = em.get_pi_metadata(observation.proposal.id)
        if program is not None:
            observation.proposal.pi_name = program['pi_name']
            observation.proposal.title = program['title']

    except Exception as e:
        logging.error('Error {} for {}'.format(e, observation.observation_id))
        tb = traceback.format_exc()
        logging.error(tb)
        return None
    logging.debug('Done update.')
    return observation


def _build_chunk_energy(chunk, filter_name, fm):

    # If n_axis=1 (as I guess it will be for all but processes GRACES
    # spectra now?) that means crpix=0.5 and the corresponding crval would
    # central_wl - bandpass/2.0 (i.e. the minimum wavelength).   It is fine
    # if you instead change crpix to 1.0.   I guess since the ‘range’ of
    # one pixel is 0.5 to 1.5.

    axis = CoordAxis1D(axis=Axis(ctype='WAVE', cunit='um'))
    ref_coord1 = RefCoord(0.5, fm.central_wl - fm.bandpass / 2.0)
    ref_coord2 = RefCoord(1.5, fm.central_wl + fm.bandpass / 2.0)
    axis.range = CoordRange1D(ref_coord1, ref_coord2)

    # DB - 14-02-19 value looks clearer (as two filters) with a space on
    # both sides of the ‘+’.
    bandpass_name = None if len(filter_name) == 0 \
        else filter_name.replace('+', ' + ')

    energy = SpectralWCS(axis=axis,
                         specsys='TOPOCENT',
                         ssyssrc='TOPOCENT',
                         ssysobs='TOPOCENT',
                         bandpass_name=bandpass_name,
                         resolving_power=fm.resolving_power)
    chunk.energy = energy
    chunk.energy_axis = 4


# values from
# https://www.gemini.edu/sciops/instruments/niri/spectroscopy/grisms
# units are? page is in microns
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


def _update_chunk_energy_niri(chunk, data_product_type, obs_id, filter_name):
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

    # https://www.gemini.edu/sciops/instruments/niri/spectroscopy/blocking-filters

    filter_md = em.get_filter_metadata(em.Inst.NIRI, filter_name)
    if filter_md is None:
        logging.warning('{}: mystery filter {} for {}'.format(
            em.Inst.NIRI, filter_name, obs_id))
        return
    filter_name = em.om.get('filter_name')
    if data_product_type == DataProductType.IMAGE:
        logging.debug('NIRI: SpectralWCS imaging for {}.'.format(obs_id))
        fm = filter_md
        fm.adjust_resolving_power()
    elif data_product_type == DataProductType.SPECTRUM:
        logging.debug('NIRI: SpectralWCS spectroscopy for {}.'.format(obs_id))
        fm = FilterMetadata('NIRI')
        fm.central_wl = filter_md.central_wl
        fm.bandpass = filter_md.bandpass
        disperser = em.om.get('disperser')
        if disperser in ['Jgrism', 'Jgrismf32', 'Hgrism', 'Hgrismf32',
                         'Kgrism', 'Kgrismf32', 'Lgrism', 'Mgrism']:
            bandpass_name = disperser[0]
            f_ratio = em.om.get('focal_plane_mask')
            logging.debug(
                'NIRI: Bandpass name is {} f_ratio is {} for {}'.format(
                    bandpass_name, f_ratio, obs_id))
            if bandpass_name in NIRI_RESOLVING_POWER:
                fm.resolving_power = \
                    NIRI_RESOLVING_POWER[bandpass_name][f_ratio]
            else:
                logging.info('NIRI: No resolving power for {}.'.format(obs_id))
        else:
            logging.info(
                'NIRI: Mystery disperser value {} for {}'.format(disperser,
                                                                 obs_id))
    else:
        raise mc.CadcException(
            'NIRI: Do not understand mode {} for {}'.format(data_product_type, obs_id))

    _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_niri')


def _update_chunk_energy_gpi(chunk, data_product_type, obs_id, filter_name):
    """GPI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_gpi')
    mc.check_param(chunk, Chunk)

    # DB - 22-02-19
    # For imaging energy WCS, use standard imaging algorithm. i.e central
    # wavelength, bandpass from SVO filter, and
    # resolution = central_wavelength/bandpass
    #
    # For energy WCS:
    #
    # naxis = extension header NAXIS2 value
    #
    # Need to hard-code the resolution from values in top table on this page:
    # https://www.gemini.edu/sciops/instruments/gpi/observing-modes-metamodes.
    # Use mean of column 4 range.  e.g. for H filter use 46.5
    # units are microns
    gpi_lookup = {'Y': (36 + 34) / 2.0,
                  'J': (39 + 35) / 2.0,
                  'H': (49 + 44) / 2.0,
                  'K1': (70 + 62) / 2.0,
                  'K2': (83 + 75) / 2.0}

    filter_md = em.get_filter_metadata(em.Inst.GPI, filter_name)
    if data_product_type == DataProductType.IMAGE:
        logging.debug('GPI: SpectralWCS imaging mode for {}.'.format(obs_id))
        fm = filter_md
    elif data_product_type == DataProductType.SPECTRUM:
        logging.debug(
            'GPI: SpectralWCS Spectroscopy mode for {}.'.format(obs_id))
        fm = FilterMetadata()
        fm.central_wl = filter_md.central_wl
        fm.bandpass = filter_md.bandpass
        if filter_name in gpi_lookup:
            fm.resolving_power = gpi_lookup[filter_name]
        else:
            raise mc.CadcException(
                'GPI: Mystery filter name {} for resolving power {}.'.format(
                    filter_name, obs_id))
    else:
        raise mc.CadcException(
            'GPI: Do not understand DataProductType {} for {}'.format(
                data_product_type, obs_id))

    _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_gpi')


def _update_chunk_energy_f2(chunk, header, data_product_type, obs_id, filter_name):
    """NIRI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_f2')
    mc.check_param(chunk, Chunk)

    filter_md = em.get_filter_metadata(em.Inst.F2, filter_name)
    if data_product_type == DataProductType.IMAGE:
        logging.debug('SpectralWCS: F2 imaging mode.')
        fm = filter_md
    elif data_product_type == DataProductType.SPECTRUM:
        logging.debug('SpectralWCS: F2 LS|Spectroscopy mode.')
        fp_mask = header.get('MASKNAME')
        mode = em.om.get('mode')
        if mode == 'LS':
            slit_width = fp_mask[0]
        else:
            # DB - 04-04-19 no way to determine slit widths used in
            # custom mask, so assume 2
            slit_width = '2'
        fm = FilterMetadata()
        fm.central_wl = filter_md.central_wl
        fm.bandpass = filter_md.bandpass
        grism_name = header.get('GRISM')
        logging.debug(
            'F2: grism name is {} fp_mask is {}'.format(grism_name, fp_mask))
        # lookup values from
        # https://www.gemini.edu/sciops/instruments/flamingos2/spectroscopy/longslit-spectroscopy
        lookup = {'1': [1300.0, 3600.0],
                  '2': [900.0, 2800.0],
                  '3': [600.0, 1600.0],
                  '4': [350.0, 1300.0],
                  '6': [130.0, 1000.0],
                  '8': [100.0, 750.0]}
        if grism_name.startswith('R3K_'):
            fm.resolving_power = lookup[slit_width][1]
        else:
            fm.resolving_power = lookup[slit_width][0]
    else:
        raise mc.CadcException(
            'Do not understand DataProductType {} for {}'.format(
                data_product_type, obs_id))

    _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_f2')


def _update_chunk_energy_general(chunk, instrument, allowable_types,
                                 data_product_type, obs_id, filter_name):
    """General chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_{}'.format(instrument))
    mc.check_param(chunk, Chunk)

    if data_product_type in allowable_types:
        logging.debug('{} Spectral WCS {} mode.'.format(instrument,
                                                        data_product_type))
        filter_md = em.get_filter_metadata(instrument, filter_name)
        _build_chunk_energy(chunk, filter_name, filter_md)
    else:
        raise mc.CadcException(
            '{} no Spectral WCS support when DataProductType {} for {}'.format(
                instrument, data_product_type, obs_id))

    logging.debug('End _update_chunk_energy_{}'.format(instrument))


# select filter_id, wavelength_central, wavelength_lower, wavelength_upper
# from gsa..gsa_filters where instrument = 'NICI'
# 0 - central
# 1 - lower
# 2 - upper
#
# dict with the barcodes stripped from the names as returned by query
# DB - 22-02-19 - units are microns
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


def _update_chunk_energy_nici(chunk, data_product_type, obs_id, filter_name):
    """NICI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_nici')
    mc.check_param(chunk, Chunk)

    filter_name = filter_name.split('_G')[0]
    filter_md = em.get_filter_metadata(em.Inst.NICI, filter_name)

    if data_product_type == DataProductType.IMAGE:
        logging.debug('NICI: SpectralWCS imaging mode for {}.'.format(obs_id))
        if filter_md is None:  # means filter_name not found
            w_max = 10.0
            w_min = 0.0
            for ii in filter_name.split('+'):
                if ii in NICI:
                    wl_max = NICI[ii][2]
                    wl_min = NICI[ii][1]
                else:
                    raise mc.CadcException(
                        'NICI: Unprepared for filter {} from {}'.format(
                            filter_name, obs_id))
                if wl_max < w_max:
                    w_max = wl_max
                if wl_min > w_min:
                    w_min = wl_min
            fm = FilterMetadata()
            fm.set_bandpass(w_max, w_min)
            fm.set_central_wl(w_max, w_min)
            fm.set_resolving_power(w_max, w_min)

        temp = em.om.get('filter_name')
        # NICI has two different bandpass names (most of the time) in its two
        # chunks.  Pat says in this case nothing will be put in the bandpass
        # name for the plane.  Add code to combine the two chunk bandpass names
        # to create a plane bandpass name only for this instrument
        _build_chunk_energy(chunk, temp, fm)
    else:
        raise mc.CadcException(
            'NICI: Do not understand DataProductType {} from {}'.format(
                data_product_type, obs_id))

    logging.debug('End _update_chunk_energy_nici')


def _update_chunk_energy_trecs(chunk, data_product_type, obs_id, filter_name):
    """TReCS-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_trecs')
    mc.check_param(chunk, Chunk)

    # Look at the json ‘disperser’ value.  If = LowRes-20" then
    # resolving power = 80.  If LowRes-10" then resolving power = 100.
    # There was a high-res mode but perhaps never used.  Again,
    # naxis = NAXIS1 header value and other wcs info as above for NIRI.  But
    # might take some string manipulation to match filter names with SVO
    # filters.

    filter_md = em.get_filter_metadata(em.Inst.TRECS, filter_name)

    if data_product_type == DataProductType.IMAGE:
        logging.debug('TRECS: imaging mode for {}.'.format(obs_id))
    elif data_product_type == DataProductType.SPECTRUM:
        logging.debug('TReCS: LS|Spectroscopy mode for {}.'.format(obs_id))
        disperser = em.om.get('disperser')
        if disperser is not None:
            if disperser == 'LowRes-20':
                filter_md.resolving_power = 80.0
            elif disperser == 'LowRes-10':
                filter_md.resolving_power = 100.0
    else:
        raise mc.CadcException(
            'Do not understand mode {} for {}'.format(data_product_type,
                                                      obs_id))

    _build_chunk_energy(chunk, filter_name, filter_md)
    logging.debug('End _update_chunk_energy_trecs')


def _update_chunk_energy_nifs(chunk, data_product_type, obs_id, filter_name):
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

    # DB - 01-03-19
    # NIFS:  this apparently uses ‘K_Short’ and ‘K_Long’ gratings (in json
    # ‘disperser’) that I think are the same as the other K grating but tuned
    # (rotated) to different wavelengths than normal.  So these should use
    # the same resolution values as the K grating lookup.

    nifs_lookup = {'Z':       {'ZJ': [1.05, 0.94, 1.15, 4990.0, 60.1]},
                   'J':       {'ZJ': [1.25, 1.15, 1.33, 6040.0, 49.6]},
                   'H':       {'JH': [1.65, 1.49, 1.80, 5290.0, 56.8]},
                   'K':       {'HK': [2.20, 1.99, 2.40, 5290.0, 56.7]},
                   'K_Short': {'HK': [2.20, 1.98, 2.41, 5290.0, 56.7]},
                   'K_Long':  {'HK': [2.20, 1.98, 2.41, 5290.0, 56.7]}}

    fm = em.get_filter_metadata(em.Inst.NIFS, filter_name)

    if data_product_type == DataProductType.SPECTRUM:
        logging.debug('NIFS: spectroscopy for {}.'.format(obs_id))
        grating = em.om.get('disperser')
        if grating in nifs_lookup:
            if filter_name in nifs_lookup[grating]:
                fm = FilterMetadata('NIFS')
                fm.set_bandpass(nifs_lookup[grating][filter_name][2],
                                nifs_lookup[grating][filter_name][1])
                fm.central_wl = em.om.get('central_wavelength')
                fm.resolving_power = nifs_lookup[grating][filter_name][3]
            else:
                raise mc.CadcException(
                    'NIFS: mystery filter_name {} for {}'.format(
                        filter_name, obs_id))
        else:
            raise mc.CadcException(
                'NIFS: mystery grating {} for {}'.format(grating, obs_id))
    elif data_product_type == DataProductType.IMAGE:
        logging.debug('NIFS: imaging for {}.'.format(obs_id))
        # DB - 01-03-19
        # NIFS images should just use the standard imaging procedure for
        # resolution (central_wavelength/bandpass).
    else:
        raise mc.CadcException(
            'NIFS: No spectroscopy information for {}'.format(obs_id))

    _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_nifs')


def _update_chunk_energy_gnirs(chunk, data_product_type, obs_id, filter_name):
    """GNIRS-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_gnirs')
    mc.check_param(chunk, Chunk)

    # DB - 07-02-19
    # Spectroscopy:
    #
    # a) long-slit
    #
    # (Note: spatial WCS info is in the primary header for GNIRS apparently)
    #
    # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy is the
    # relevant web page.  Grating, Order-Band, Blocking Filter Range and the
    # two Resolving Power columns (one for each camera) are the important
    # columns.
    #
    # Need to know:
    #
    # grating: from json ‘disperser’ (need 10, 32 or 111 numbers only)
    # Note: if disperser contains string “XD” then the spectrum is cross
    # dispersed.  See b) below.
    #
    # filter: from json ‘filter_name’
    # camera: Short or Long substring from json ‘camera’ value (i.e. “Blue”
    # or “Red” aren’t important)
    # central wavelength: from json ‘central_wavelength’
    # NAXIS2 extension header value
    #
    # crval = central_wavelength
    #
    # use ‘Blocking Filter Range’ for the appropriate filter
    # (e.g. X, J, H...) to determine min/max wavelengths as the rough
    # estimate of the wavelength coverage.

    # b) cross-dispersed (when ‘disperser’ contains ‘?XD’ string)
    #
    # Wavelength ranges given here:
    # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/crossdispersed-xd-spectroscopy/xd-prisms
    #
    # You need the ‘Short’ or ‘Long’ from the json ‘camera’ value and the
    # ‘SXD’ or ‘LXD’ from the ‘disperser’ value to look up the coverage.
    #
    # Imaging:
    #
    # Not done very often.  Filter info is here:
    # https://www.gemini.edu/sciops/instruments/gnirs/imaging.  For filter
    # ID’s I think json ‘filter_name’ value of J_(order_5) corresponds to
    # “J (order blocking)” in this table.  ‘filter_name’ value of J
    # corresponds to “J (Mauna Kea)“.

    # DB - 01-03-19
    # Resolution is in the last two columns of the table here, a different
    # value for each grating/filter combination as well as camera value:
    # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy.
    # Add that resolution for each grating/filter (and camera) combination
    # sort of as \for the bandpasses.  Already know the grating value.
    # Then, for both long-slit and cross-dispersed spectroscopy need to look
    # at json ‘camera’ (or header CAMERA) value (does it contain ‘short’ or
    # ‘long’) and the json focal_plane_mask value (or SLIT header value;
    # you need only the numeric value at the start of the string).  The
    # camera value tells you which column to look under for the resolution.
    # The best estimate of the resolution, R, is then given by:
    #
    # R = 0.3 x tabulated value / slit width  (for ‘short’ camera observations)
    # R = 0.1 x tabulated value / slit width  (for ‘long’ camera observations)
    #
    # i.e. the resolution gets lower when a wider slit is in place
    # (the tabulated values are for 0.3" and 0.1" slits for the short/long
    # cameras respectively, hence the numbers in the start of each equation).
    # Basically the camera is producing an image of the slit at each
    # wavelength but if the slit is wider the dispersed image is also wider
    # and so different colours are blended together more so the resolution
    # gets worse as you open up the slit (but you get more light through a
    # wider slit so it’s a compromise between throughput and spectral
    # resolution).

    fm = FilterMetadata('GNIRS')
    if data_product_type == DataProductType.SPECTRUM:
        logging.debug(
            'gnirs: SpectralWCS Spectroscopy mode for {}.'.format(obs_id))
        fm.central_wl = em.om.get('central_wavelength')
        disperser = em.om.get('disperser')
        if 'XD' in disperser:
            logging.error('gnirs: cross dispersed mode.')
            # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/
            # crossdispersed-xd-spectroscopy/xd-prisms
            # 0 = lower
            # 1 = upper
            # 2 = spectral resolution with 2-pix wide slit
            xd_mode = {'SB+SXD': [0.9, 2.5, 2.0],
                       'LB+LXD': [0.9, 2.5, 2.3],
                       'LB+SXD': [1.2, 2.5, 2.5]}
            lookup = None
            coverage = disperser[-3:]
            camera = em.om.get('camera')
            if camera.startswith('Short'):
                lookup = '{}+{}'.format('SB', coverage)
            elif camera.startswith('Long'):
                lookup = '{}+{}'.format('LB', coverage)
            if lookup is None:
                raise mc.CadcException(
                    'gnirs: Do not understand xd mode {} {}'
                        .format(camera, coverage))
            bounds = xd_mode[lookup]

            focal_plane_mask = em.om.get('focal_plane_mask')
            slit_width = 1.0
            if 'arcsec' in focal_plane_mask:
                slit_width = float(focal_plane_mask.split('arcsec')[0])
            if camera.startswith('Long'):
                fm.resolving_power = _calc_gnirs_resolving_power(
                    0.1, xd_mode[lookup][2], slit_width)
            elif camera.startswith('Short'):
                fm.resolving_power = _calc_gnirs_resolving_power(
                    0.3, xd_mode[lookup][2], slit_width)

        else:
            logging.debug('gnirs: long slit mode for {}.'.format(obs_id))
            # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy
            # 0 - min wavelength
            # 1 - max wavelength
            # 2 - 'short' camera resolution
            # 3 - 'long' camera resolution
            # 4 - Since November 2012 and for the cross-dispersed mode with
            # the 2 pix wide slit only resolving powers are somewhat lower, as
            # follows: X-1400; J-1400, H-1400; K-1300
            long_slit_mode = {'11': {'X': [1.03, 1.17, 570, 2100],
                                     'J': [1.17, 1.37, 570, 1600],
                                     'H': [1.47, 1.80, 570, 1700],
                                     'K': [1.91, 2.49, 570, 1700],
                                     'L': [2.80, 4.20, 570, 1800],
                                     'M': [4.40, 6.00, 570, 1200]},
                              '32': {'X': [1.03, 1.17, 1700, 5100, 1400],
                                     'J': [1.17, 1.37, 1600, 4800, 1400],
                                     'H': [1.49, 1.80, 1700, 5100, 1400],
                                     'K': [1.91, 2.49, 1700, 5100, 1300],
                                     'L': [2.80, 4.20, 1800, 5400, 1800],
                                     'M': [4.40, 6.00, 1240, 3700, 1240]},
                              '111': {'X': [1.03, 1.17, 6600, 17800],
                                      'J': [1.17, 1.37, 7200, 17000],
                                      'H': [1.49, 1.80, 5900, 17800],
                                      'K': [1.91, 2.49, 5900, 17800],
                                      'L': [2.80, 4.20, 6400, 19000],
                                      'M': [4.40, 6.00, 4300, 12800]}}
            bandpass = filter_name[0]
            grating = disperser.split('_')[0]
            bounds = long_slit_mode[grating][bandpass]
            if not bounds:
                raise mc.CadcException(
                    'gnirs: Do not understand long slit mode with {} {}'
                        .format(bandpass, bounds))
            camera = em.om.get('camera')
            focal_plane_mask = em.om.get('focal_plane_mask')
            slit_width = 1.0
            if 'arcsec' in focal_plane_mask:
                slit_width = float(focal_plane_mask.split('arcsec')[0])
            if camera.startswith('Long'):
                fm.resolving_power = _calc_gnirs_resolving_power(
                    0.1, long_slit_mode[grating][bandpass][2], slit_width)
            elif camera.startswith('Short'):
                date_time = ac.get_datetime(em.om.get('ut_datetime'))
                if date_time > ac.get_datetime('2012-11-01T00:00:00'):
                    fm.resolving_power = _calc_gnirs_resolving_power(
                        0.3, long_slit_mode[grating][bandpass][4], slit_width)
                else:
                    fm.resolving_power = _calc_gnirs_resolving_power(
                        0.3, long_slit_mode[grating][bandpass][3], slit_width)
            else:
                raise mc.CadcException(
                    'GNIRS: Mystery camera value {} for {}'.format(camera,
                                                                   obs_id))

        fm.set_bandpass(bounds[1], bounds[0])

    elif data_product_type == DataProductType.IMAGE:
        logging.debug('gnirs: SpectralWCS imaging mode for {}.'.format(obs_id))
        # https://www.gemini.edu/sciops/instruments/gnirs/imaging
        # 0 - min wavelength
        # 1 - max wavelength
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

        fm.set_central_wl(bounds[1], bounds[0])
        fm.set_bandpass(bounds[1], bounds[0])
    else:
        raise mc.CadcException(
            'GNIRS: Unexpected DataProductType {} for {}'.format(
                data_product_type, obs_id))
    _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_gnirs')


def _calc_gnirs_resolving_power(ratio, initial_value, slit_width):
    return ratio * initial_value / slit_width


# select filter_id, wavelength_central, wavelength_lower, wavelength_upper
# from gsa..gsa_filters where instrument = 'NICI'
# 0 - central
# 1 - lower
# 2 - upper
#
# units are microns
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


def _update_chunk_energy_phoenix(chunk, data_product_type, obs_id, filter_name):
    """Phoenix-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_phoenix')
    mc.check_param(chunk, Chunk)

    # DB - 12-02-19 - Phoenix should be the same as TReCS but uses NAXIS2 for
    # the length of the dispersion axis.

    # DB - 12-02-19 - Note that the parenthetical numbers
    # after the Phoenix filter names (in the header) indicates the
    # filter wheel slot the filter is in and may occasionally change
    # so should be disregarded.
    if len(filter_name) > 0:
        filter_name = filter_name.split()[0]

    logging.debug(
        'Phoenix: filter_name is {} for {}'.format(filter_name, obs_id))

    fm = FilterMetadata('Phoenix')
    if data_product_type in [DataProductType.SPECTRUM,
                             DataProductType.IMAGE]:
        logging.debug(
            'Phoenix: DataProductType {} for {}.'.format(data_product_type,
                                                         obs_id))
        if filter_name in PHOENIX:
            fm.set_bandpass(PHOENIX[filter_name][2], PHOENIX[filter_name][1])
            fm.central_wl = PHOENIX[filter_name][0]
        elif len(filter_name) == 0:
            fm.set_bandpass(10.0, 0.0)
            fm.set_central_wl(10.0, 0.0)
        else:
            raise mc.CadcException(
                'Phoenix: mystery filter name {} for {}'.format(
                    filter_name, obs_id))
    else:
        raise mc.CadcException(
            'Phoenix: Unsupported DataProductType {} for {}'.format(
                data_product_type, obs_id))

    _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_phoenix')


def _update_chunk_energy_flamingos(chunk, data_product_type, obs_id, filter_name):
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

    # spectral wcs units are microns, values from Table 7 are angstroms.
    # The conversion is here.

    # DB - 21-02-19 - JH central=1.45 um, FWHM=0.95 um
    # 0 = central wavelength
    # 1 = FWHM
    lookup = {'JH': [1.45, 0.95],
              'HK': [(2.7588 + 1.0347) / 2.0, (2.7588 - 1.0347)]}

    fm = FilterMetadata()
    if filter_name in lookup:
        fm.central_wl = lookup[filter_name][0]
        fm.bandpass = lookup[filter_name][1]
    else:
        fm = em.get_filter_metadata(em.Inst.FLAMINGOS, filter_name)

    if data_product_type == DataProductType.SPECTRUM:
        logging.debug('Flamingos: SpectralWCS for {}.'.format(obs_id))
        fm.resolving_power = 1300.0
    elif data_product_type == DataProductType.IMAGE:
        logging.debug(
            'Flamingos: SpectralWCS imaging mode for {}.'.format(obs_id))
    else:
        raise mc.CadcException(
            'Flamingos: mystery data product type {} for {}'.format(
                data_product_type, obs_id))

    _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_flamingos')


def _update_chunk_energy_hokupaa(chunk, data_product_type, obs_id, filter_name):
    """hokupaa-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_hokupaa')
    mc.check_param(chunk, Chunk)

    # DB - 18-01-19 - Note: it is always imaging for Hokupa'a + QUIRC so that
    # could be hardcoded.

    # DB - 20-02-19 - This page has a table of QUIRC (the camera part of
    # Hokupa’a + QUIRC) filter names and bandpasses (in microns).
    # Use this as a lookup for central wavelengths and bandpass
    # http://www.gemini.edu/sciops/instruments/uhaos/uhaosQuirc.html

    # J+CO is "Dark Position", units are microns
    # 0 - central wavelength
    # 1 - bandpass
    hokupaa_lookup = {'J': [1.25, 0.171],
                      'H': [1.65, 0.296],
                      'K': [2.2, 0.336],
                      'K\'': [2.12, 0.41],
                      'H+K notch': [1.8, 0.7],
                      'methane low': [1.56, 120.0],
                      'methane high': [1.71, 120.0],
                      'FeII': [1.65, 170.0],
                      'HeI': [2.06, 30.0],
                      '1-0 S(1) H2': [2.12, 23.0],
                      'H Br(gamma)': [2.166, 150.0],
                      'K-continuum': [2.26, 60.0],
                      'CO': [2.29, 20.0]
                      }

    if filter_name not in hokupaa_lookup:
        raise mc.CadcException(
            'hokupaa: Mystery filter {} for {}'.format(filter_name, obs_id))
    if data_product_type == DataProductType.IMAGE:
        logging.debug(
            'hokupaa: SpectralWCS imaging mode for {}.'.format(obs_id))
        fm = FilterMetadata()
        fm.central_wl = hokupaa_lookup[filter_name][0]
        fm.bandpass = hokupaa_lookup[filter_name][1]
    else:
        raise mc.CadcException(
            'hokupaa: mystery data product type {} for {}'.format(
                data_product_type, obs_id))
    _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_hokupaa')


def _update_chunk_energy_oscir(chunk, data_product_type, obs_id, filter_name):
    """oscir-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_oscir')
    mc.check_param(chunk, Chunk)

    # Filter info here:
    # http://www.gemini.edu/sciops/instruments/oscir/oscirFilterList.html
    # No filter provided in json; use FILTER keyword.
    # e.g. ‘S_12.5 (-11775)’ = ‘12.5’ in table.
    # It looks like only the 'r' files have filter ids.

    # 0 - central wavelenth
    # 1 - bandpass
    # units are microns
    oscir_lookup = {'7.9': [7.91, 0.755],
                    '8.8': [8.81, 0.871],
                    '9.8': [9.80, 0.952],
                    '10.3': [10.27, 1.103],
                    '11.7': [11.70, 1.110],
                    '12.5': [12.49, 1.156],
                    'N': [10.75, 5.230],
                    'IHW': [18.17, 1.651],
                    'Q3': [20.8, 1.650]}

    temp = filter_name
    filter_name = temp.split('_')[0]
    if filter_name not in oscir_lookup:
        raise mc.CadcException(
            'oscir: Mystery FILTER keyword {} for {}'.format(
                filter_name, obs_id))
    if data_product_type == DataProductType.IMAGE:
        logging.debug(
            'oscir: SpectralWCS imaging mode for {}.'.format(obs_id))
        fm = FilterMetadata()
        fm.central_wl = oscir_lookup[filter_name][0]
        fm.bandpass = oscir_lookup[filter_name][1]
        _build_chunk_energy(chunk, filter_name, fm)
    else:
        raise mc.CadcException(
            'oscir: mystery data product type {} for {}'.format(
                data_product_type, obs_id))

    logging.debug('End _update_chunk_energy_oscir')


def _update_chunk_energy_bhros(chunk, header, data_product_type, obs_id):
    """bhros-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_bhros')
    mc.check_param(chunk, Chunk)
    # DB - 20-02-19 - There were bHROS filters but I don’t think they were
    # used during the very limited lifetime of the instrument.  No info
    # in the headers either.
    #
    # bHROS spectral resolution should be approximately 150,000/x-binning
    # value. json returns a “detector_binning”: “1x1" value where the 1x1
    # indicates no binning in x or y (for this example).  Could be 2x1,
    # 4x2, etc.  Binning is determined from header keyword:
    # CCDSUM  = ‘1 1     ’           / CCD pixel summing
    #
    # The approximate central wavelength is json ‘central_wavelength’ value.
    # Unfortunately the wavelength coverage is not straightforward.
    # See http://www.gemini.edu/sciops/instruments/hros/hrosDispersion.html.
    # The CCD did not cover the entire spectrum so only a subset of the
    # entire optical spectral region was observed.

    # Use central wavelength in microns and +/- 0.2 microns as a better
    # guess-timate rather than imply that entire spectrum is present.

    if data_product_type == DataProductType.SPECTRUM:
        logging.debug('bhros: SpectralWCS spectroscopy for {}.'.format(obs_id))
        fm = FilterMetadata()
        fm.central_wl = em.om.get('central_wavelength')
        fm.adjust_bandpass(0.2)
        fm.resolving_power = 150000.0
        ccd_sum = header.get('CCDSUM')
        if ccd_sum is not None:
            temp = float(ccd_sum.split()[1])
            fm.resolving_power = 150000.0 / temp
        _build_chunk_energy(chunk, '', fm)
    else:
        raise mc.CadcException(
            'bhros: mystery data product type {} for {}'.format(
                data_product_type, obs_id))

    logging.debug('End _update_chunk_energy_bhros')


def _update_chunk_energy_graces(chunk, data_product_type, obs_id):
    """graces-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_graces')
    mc.check_param(chunk, Chunk)

    # DB - 22-02-19  Axis 2 is the spectral axis.   Resolution is
    # 67,500 divided by second number in json ‘detector_binning’
    # value of ‘N x N’.  The central wavelength in json
    # is always set to 0.7 microns and range is from about 0.4 to
    # 1.0 microns.  So use 0.7 as the crval, crpix = naxis2/2.0 and
    # delta = (1.0 - 0.4)/naxis2.  Just a kludge to get spectral
    # range correct for an echelle spectrograph.

    if data_product_type == DataProductType.SPECTRUM:
        logging.debug('graces: SpectralWCS spectroscopy for {}.'.format(obs_id))
        fm = FilterMetadata()
        fm.central_wl = em.om.get('central_wavelength')
        fm.set_bandpass(1.0, 0.4)
        ccd_sum = em.om.get('detector_binning')
        fm.resolving_power = 67500.0
        if ccd_sum is not None:
            temp = float(ccd_sum.split('x')[1])
            fm.resolving_power = 67500.0 / temp
        _build_chunk_energy(chunk, '', fm)
    else:
        raise mc.CadcException(
            'graces: mystery data product type {} for {}'.format(
                data_product_type, obs_id))

    logging.debug('End _update_chunk_energy_graces')


def _update_chunk_energy_gmos(chunk, data_product_type, obs_id, filter_name,
                              instrument):
    logging.debug('Begin _update_chunk_energy_gmos')

    GMOS_RESOLVING_POWER = {
        'B1200': 3744.0,
        'R831': 4396.0,
        'B600': 1688.0,
        'R600': 3744.0,
        'R400': 1918.0,
        'R150': 631.0
    }
    # 0 = min
    # 1 = max
    # units are microns
    lookup = {'GG455': [0.46000, 1.10000],
              'OG515': [0.52000, 1.10000],
              'RG610': [0.61500, 1.10000],
              'RG780': [0.07800, 1.10000],
              'open': [0.35000, 1.10000]}

    filter_md = em.get_filter_metadata(instrument, filter_name)
    if filter_md is None:  # means filter_name not found
        w_max = 10.0
        w_min = 0.0
        for ii in filter_name.split('+'):
            if 'Hartmann' in ii:
                continue
            elif ii in lookup:
                wl_max = lookup[ii][1]
                wl_min = lookup[ii][0]
            else:
                raise mc.CadcException(
                    '{}: Mystery filter {} from {}'.format(
                        instrument, filter_name, obs_id))
            if wl_max < w_max:
                w_max = wl_max
            if wl_min > w_min:
                w_min = wl_min
        filter_md = FilterMetadata()
        filter_md.set_bandpass(w_max, w_min)
        filter_md.set_central_wl(w_max, w_min)
        filter_md.set_resolving_power(w_max, w_min)

    if data_product_type == DataProductType.SPECTRUM:
        logging.debug('gmos: SpectralWCS spectroscopy for {}.'.format(obs_id))
        if math.isclose(filter_md.central_wl, 0.0):
            logging.info(
                'gmos: no spectral wcs, central wavelength is {} for {}'.format(
                    filter_md.central_wl, obs_id))
            return
        disperser = em.om.get('disperser')
        fm = FilterMetadata()
        fm.central_wl = filter_md.central_wl
        fm.bandpass = filter_md.bandpass
        if disperser in GMOS_RESOLVING_POWER:
            fm.resolving_power = GMOS_RESOLVING_POWER[disperser]
        else:
            raise mc.CadcException(
                'gmos: mystery disperser {} for {}'.format(disperser, obs_id))
    elif data_product_type == DataProductType.IMAGE:
        logging.debug('gmos: SpectralWCS imaging for {}.'.format(obs_id))
        fm = filter_md
    else:
        raise mc.CadcException(
            'gmos: mystery data product type {} for {}'.format(
                data_product_type, obs_id))
    _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_gmos')


def _reset_energy(observation_type, data_label, instrument):
    """
    Return True if there should be no energy WCS information created at
    the chunk level.

    :param observation_type from the parent Observation instance.
    :param data_label str for useful logging information only.
    """
    result = False
    om_filter_name = em.om.get('filter_name')

    if ((observation_type is not None and ((observation_type == 'DARK') or
         (instrument in [em.Inst.GMOS, em.Inst.GMOSN, em.Inst.GMOSS] and
          observation_type == 'BIAS'))) or
        (om_filter_name is not None and ('blank' in om_filter_name or
                                         'Blank' in om_filter_name))):
        logging.info(
            'No chunk energy for {} obs type {} filter name {}'.format(
                data_label, observation_type, om_filter_name))
        result = True
    return result


def get_filter_name(primary_header, header, obs_id, instrument=None):
    """
    Create the filter names.

    :param primary_header: The zero-th FITS header.
    :param header: The FITS header for the current extension.
    :param obs_id: The observation for which filter name search is happening
        (used for logging only)
    :param instrument: For instrument-specific behaviour.
    :return: The filter names, or None if none found.
    """
    if instrument == em.Inst.GRACES:
        # filter name not used in spectral WCS calculation
        return None

    header_filters = []

    filter_name = None
    #
    # NIRI - prefer header keywords
    #
    # NICI - use filter names from headers, because there's a different
    # filter/header, and the JSON summary value obfuscates that
    #
    if instrument not in [em.Inst.NIRI, em.Inst.NICI]:
        filter_name = em.om.get('filter_name')
    if (filter_name is None or
            (filter_name is not None and len(filter_name.strip()) == 0)):
        # DB - 04-02-19 - strip out anything with 'pupil' as it doesn't affect
        # energy transmission
        filters2ignore = ['open', 'invalid', 'pupil']
        lookup = 'FILTER'
        if instrument == em.Inst.PHOENIX:
            lookup = 'FILT_POS'
            search_header = header
        elif instrument == em.Inst.NICI:
            search_header = header
        else:
            search_header = primary_header
        for key in search_header.keys():
            if lookup in key:
                value = search_header.get(key).lower()
                ignore = False
                for ii in filters2ignore:
                    if ii.startswith(value) or value.startswith(ii):
                        ignore = True
                        break
                if ignore:
                    continue
                else:
                    header_filters.append(search_header.get(key).strip())
            filter_name = '+'.join(header_filters)
    logging.info(
        'Filter names are {} for instrument {} in {}'.format(filter_name,
                                                             instrument,
                                                             obs_id))
    return filter_name


def _update_position_from_zeroth_header(artifact, headers, log_id):
    """Make the 0th header spatial WCS the WCS for all the
    chunks."""
    primary_header = headers[0]

    # naxis values are only available from extensions

    primary_header['NAXIS1'] = headers[-1].get('NAXIS1')
    primary_header['NAXIS2'] = headers[-1].get('NAXIS2')

    wcs_parser = WcsParser(primary_header, log_id, 0)
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


def _update_chunk_position_flamingos(chunk, header, obs_id):
    # DB - I see nothing in astropy that will do a transformation from crota
    # form to CD matrix, but this is it:

    # cd1_1 = cdelt1 * cos (crota1)
    # cd1_2 = -cdelt2 * sin (crota1)
    # cd2_1 = cdelt1 * sin (crota1)
    # cd2_2 = cdelt2 * cos (crota1)

    # Note that there is not a crota2 keyword (it would have the same value
    # as crota1 if it existed)
    if (chunk is not None and chunk.position is not None and
            chunk.position.axis is not None and
            chunk.position.axis.function is not None):
        c_delt1 = header.get('CDELT1')
        c_delt2 = header.get('CDELT2')
        c_rota1 = header.get('CROTA1')
        if c_delt1 is not None and c_delt2 is not None and c_rota1 is not None:
            chunk.position.axis.function.cd11 = c_delt1 * math.cos(c_rota1)
            chunk.position.axis.function.cd12 = -c_delt2 * math.sin(c_rota1)
            chunk.position.axis.function.cd21 = c_delt1 * math.sin(c_rota1)
            chunk.position.axis.function.cd22 = c_delt2 * math.cos(c_rota1)
        else:
            logging.info(
                'FLAMINGOS: Missing spatial wcs inputs for {}'.format(obs_id))
            chunk.position.axis.function.cd11 = None
            chunk.position.axis.function.cd12 = None
            chunk.position.axis.function.cd21 = None
            chunk.position.axis.function.cd22 = None
    else:
        logging.info('FLAMINGOS: Missing spatial wcs for {}'.format(obs_id))


def _update_chunk_position(chunk, header, instrument, extension, obs_id):
    logging.debug('Begin _update_chunk_position')
    mc.check_param(chunk, Chunk)

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
    #
    # DB - Hokupa'a+QUIRC
    # This is an imager, so not a single big pixel.  So build a CD matrix
    # (ignoring any possible camera rotation to start - not sure it can
    # be determined).  Use RA/Dec in header for CRVALs (no json values
    # returned I think) but you have to convert these from HH:MM:SS etc.
    # format to degrees.  CRPIX1/2 = NAXIS1/2 divided by 2.  PIXSCALE in
    # header (both primary and extension) give plate scale that is
    # fixed at 0.01998 arcsec/pixel.
    # CD1_1 = CD2_2 = PIXSCALE/3600.0 (to convert to degrees).
    # CD1_2 = CD2_1 = 0.0.
    #
    # DB - 20-02-19 - OSCIR
    # NAXIS1/2 give number of pixels so CRPIX1/2 = NAXIS1/2 divided by 2.
    # json RA/Dec are bogus.  Need to use RA_TEL and DEC_TEL and convert
    # to degrees and use these for CRVAL1/2 values.  (RA_BASE and DEC_BASE
    # values in degrees don’t quite agree with RA_TEL and DEC_TEL...)
    # Use PIXSCALE= ‘0.089’ value to build CD matrix.
    # So same as for Hokupa’a:  CD1_1 = CD2_2 = PIXSCALE/3600.0.
    # CD1_2 = CD2_1 = 0.0
    #
    # DB - 22-02-19 - GPI
    # WCS info is garbage in header, so for spatial WCS (for both image and
    # spectrum):
    #
    # crval1 = json ra value
    # crval2 = json dec value
    # naxis1 = extension header NAXIS1 value
    # naxis2 = extension header NAXIS2 value
    # crpix1 = naxis1/2.0
    # crpix2 = naxis2/2.0
    #
    # FOV is 2.8" x 2.8" on each side or 2.8/3600.0 degrees
    # cd1_1 = 2.8/(3600.0 * naxis1)
    # cd2_2 = 2.8/(3600.0 * naxis2)
    # cd1_2 = cd2_1 = 0.0 since we don’t know rotation value

    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CRVAL1'] = get_ra(header)
    header['CRVAL2'] = get_dec(header)
    if instrument != em.Inst.GPI:
        header['CDELT1'] = RADIUS_LOOKUP[instrument]
        header['CDELT2'] = RADIUS_LOOKUP[instrument]
        header['CROTA1'] = 0.0
    if instrument not in [em.Inst.OSCIR, em.Inst.GPI]:
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
    header['CRPIX1'] = get_crpix1(header)
    header['CRPIX2'] = get_crpix2(header)
    header['CD1_1'] = get_cd11(header)
    header['CD1_2'] = 0.0
    header['CD2_1'] = 0.0
    header['CD2_2'] = get_cd22(header)

    wcs_parser = WcsParser(header, obs_id, extension)
    if chunk is None:
        chunk = Chunk()
    wcs_parser.augment_position(chunk)
    chunk.position_axis_1 = 1
    chunk.position_axis_2 = 2

    if instrument == em.Inst.OSCIR:
        chunk.position.coordsys = header.get('FRAMEPA')
        chunk.position.equinox = float(header.get('EQUINOX'))
    elif instrument == em.Inst.BHROS:
        chunk.position.coordsys = header.get('TRKFRAME')
        chunk.position.equinox = float(header.get('TRKEQUIN'))
    elif instrument == em.Inst.GPI:
        chunk.position.coordsys = header.get('RADESYS')

    logging.debug('End _update_chunk_position')


def _update_chunk_time_f2(chunk, obs_id):
    """F2 FITS files have a CD3_3 element that's not supported by fits2caom2,
    so using the blueprint will not work to adjust that value. Set delta
    specifically here."""
    logging.debug('Begin _update_chunk_time_f2 {}'.format(obs_id))
    mc.check_param(chunk, Chunk)
    if (chunk.time is not None and chunk.time.axis is not None and
            chunk.time.axis.function is not None):
        chunk.time.axis.function.delta = get_time_delta(None)
        logging.info('F2: Updated time delta for {}'.format(obs_id))
    logging.debug('End _update_chunk_time_f2 {}'.format(obs_id))


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
