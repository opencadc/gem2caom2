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
import traceback

from astropy import units
from astropy.coordinates import SkyCoord

from caom2 import Observation, ObservationIntentType, DataProductType
from caom2 import CalibrationLevel, TargetType, ProductType, Chunk, Axis
from caom2 import SpectralWCS, CoordAxis1D, CoordFunction1D, RefCoord
from caom2 import TypedList, CoordRange1D, CompositeObservation, Algorithm
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2utils import WcsParser
from caom2pipe import manage_composable as mc
from caom2pipe import execute_composable as ec
from caom2pipe import caom_composable as cc
from caom2pipe import astro_composable as ac

import gem2caom2.external_metadata as em
from gem2caom2.gem_name import GemName, COLLECTION, ARCHIVE, SCHEME
from gem2caom2.svofps import FilterMetadata
from gem2caom2.obs_file_relationship import GemObsFileRelationship

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
# NIFS - DB - 04-03-19 - hard-code 3" FOV
# TEXES - DB - 07-03-19 - cd11=cd22 = 5.0/3600.0
RADIUS_LOOKUP = {em.Inst.GPI: 2.8 / 3600.0,  # units are arcseconds
                 em.Inst.GRACES: 1.2 / 3600.0,
                 em.Inst.PHOENIX: 5.0 / 3600.0,
                 em.Inst.OSCIR: 0.0890 / 3600.0,
                 em.Inst.HOKUPAA: 4.0 / 3600.0,
                 em.Inst.BHROS: 0.9 / 3600.0,
                 em.Inst.NIFS: 3.0 / 3600.0,
                 em.Inst.TEXES: 5.0 / 3600.0}


def get_time_delta(header):
    """
    Calculate the Time WCS delta.

    :param header: The FITS header for the current extension.
    :return: The Time delta, or None if none found.
    """
    exptime = get_exposure(header)
    if _is_gmos_mask(header):
        # DB - 05-03-19 - delta hardcoded to 0
        exptime = 0.0
    if exptime is None:
        return None
    return float(exptime) / (24.0 * 3600.0)


def get_calibration_level(uri):
    result = CalibrationLevel.RAW_STANDARD
    reduction = em.om.get('reduction')
    instrument = _get_instrument()
    if ((reduction is not None and
         (('PROCESSED' in reduction) or ('PREPARED' in reduction))) or
            (instrument is em.Inst.TEXES and
             ('_red' in uri.lower() or '_sum' in uri.lower())) or
            (instrument is em.Inst.PHOENIX and
             ec.CaomName(uri.lower()).file_id.startswith('p'))):
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

    logging.debug(
        'obs type is {} obs class is {} for '.format(obs_type, obs_class,
                                                     _get_data_label()))
    if obs_type is not None and obs_type == 'MASK':
        result = ProductType.AUXILIARY
    elif obs_class is None:
        # 'unknown' is the value of observation_class for CIRPASS
        if obs_type is not None and obs_type in ['OBJECT', 'unknown']:
            obs_id = header.get('DATALAB')
            if obs_id is not None and 'CAL' in obs_id:
                result = ProductType.CALIBRATION
            else:
                result = ProductType.SCIENCE
        else:
            instrument = _get_instrument()
            if instrument is not None:
                if instrument is em.Inst.PHOENIX:
                    if _is_phoenix_calibration(header):
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                elif instrument is em.Inst.OSCIR:
                    if _is_oscir_calibration(header):
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                elif instrument is em.Inst.FLAMINGOS:
                    if _is_flamingos_calibration():
                        result = ProductType.CALIBRATION
                    else:
                        result = ProductType.SCIENCE
                elif instrument is em.Inst.HOKUPAA:
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
    if instrument is em.Inst.HOKUPAA:
        result = _get_pix_scale(header)
    elif instrument in [em.Inst.OSCIR, em.Inst.TEXES]:
        result = RADIUS_LOOKUP[instrument]
    elif instrument in [em.Inst.GPI, em.Inst.NIFS]:
        # DB - 05-03-19 - NIFS needs a division by NAXIS1/2 for the
        # cdelta1/2 calculations.
        result = RADIUS_LOOKUP[instrument]/header.get('NAXIS1')
    elif instrument is em.Inst.CIRPASS:
        # DB - 06-03-19
        # FOV is fixed at two possible values and has no bearing on NAXIS1/2
        # values. See
        # http://www.gemini.edu/sciops/instruments/cirpass/cirpassIFU.html.
        # 499 lenslets cover the FOV: about 33 along one axis and 15 along the
        # other.
        #
        # LENS_SCL determines the scale/lenslet:  0.36 or 0.25 (arcseconds per lens)
        #     if 0.36 then FOV is 13.0" x 4.7" (RA and Dec)
        #     if 0.25 then FOV is 9.3" x 3.5"
        #     cd11 = LENS_SCL/3600.0
        #     cd22 = LENS_SCL/3600.0
        lens_scl = header.get('LENS_SCL')
        result = float(lens_scl) / 3600.0
    else:
        cdelt1 = header.get('CDELT1')
        if cdelt1 is None:
            result = RADIUS_LOOKUP[instrument]
        else:
            result = cdelt1
    return result


def get_cd22(header):
    instrument = _get_instrument()
    if instrument is em.Inst.HOKUPAA:
        result = _get_pix_scale(header)
    elif instrument in [em.Inst.OSCIR, em.Inst.TEXES]:
        result = RADIUS_LOOKUP[instrument]
    elif instrument in [em.Inst.GPI, em.Inst.NIFS]:
        # DB - 05-03-19 - NIFS needs a division by NAXIS1/2 for the
        # cdelta1/2 calculations.
        result = RADIUS_LOOKUP[instrument]/header.get('NAXIS2')
    elif instrument is em.Inst.CIRPASS:
        lens_scl = header.get('LENS_SCL')
        result = float(lens_scl) / 3600.0
    else:
        cdelt2 = header.get('CDELT2')
        if cdelt2 is None:
            result = RADIUS_LOOKUP[instrument]
        else:
            result = cdelt2
    return result


def get_crpix1(header):
    instrument = _get_instrument()
    if instrument is em.Inst.HOKUPAA:
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
    if instrument is em.Inst.HOKUPAA:
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
    if instrument is em.Inst.FLAMINGOS:
        result, ignore = _get_flamingos_mode(header)
    elif instrument is em.Inst.GNIRS:
        # DB - 03-04-19
        # if disperser for GNIRS = MIRROR
        # then it’s an image, otherwise a spectrum.
        disperser = em.om.get('disperser')
        if disperser == 'MIRROR':
            result = DataProductType.IMAGE
        else:
            result = DataProductType.SPECTRUM
    elif instrument in [em.Inst.CIRPASS, em.Inst.TEXES]:
        result = DataProductType.SPECTRUM
    else:
        mode = em.om.get('mode')
        obs_type = _get_obs_type(header)
        if mode is None:
            raise mc.CadcException('No mode information found for {}'.format(
                em.om.get('filename')))
        elif instrument is em.Inst.GPI:
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
    if instrument is em.Inst.HOKUPAA:
        ra, dec = _get_sky_coord(header, 'RA', 'DEC')
        result = dec
    elif instrument is em.Inst.OSCIR:
        ra, dec = _get_sky_coord(header, 'RA_TEL', 'DEC_TEL')
        result = dec
    elif instrument in [em.Inst.BHROS, em.Inst.NIFS, em.Inst.TEXES]:
        # bHROS, TEXES ra/dec not in json
        # NIFS ra/dec not reliable in json
        result = header.get('DEC')
    elif instrument is em.Inst.CIRPASS:
        # DB - 06-03-19 - Must use FITS header info for most WCS info
        ra, dec = _get_sky_coord(header, 'TEL_RA', 'TEL_DEC')
        result = dec
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
    if instrument is em.Inst.OSCIR:
        # DB - 20-02-19 - json ‘exposure_time’ is in minutes, so multiply
        # by 60.0.
        result = result * 60.0
    elif instrument in [em.Inst.FLAMINGOS, em.Inst.CIRPASS]:
        # exposure_time is null in the JSON
        # DB - 06-03-19 Use value of header keyword EXP_TIME for exposure time.
        result = mc.to_float(header.get('EXP_TIME'))
    elif instrument is em.Inst.TEXES:
        # DB - 07-03-19 -  Use header ‘OBSTIME’ value as exposure time in
        # seconds.
        result = mc.to_float(header.get('OBSTIME'))
    elif _is_gmos_mask(header):
        result = 0.0
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
                  'flat', 'ARC', 'Domeflat', 'DARK', 'dark', 'gcal']
    dl = header.get('DATALAB')
    lookup = _get_obs_class(header)
    logging.debug('observation_class is {} for {}'.format(lookup, dl))
    if lookup is None:
        type_lookup = _get_obs_type(header)
        logging.debug('observation_type is {} for {}'.format(type_lookup, dl))
        if type_lookup is None:
            data_label = _get_data_label()
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
                instrument = _get_instrument()
                if instrument is not None:
                    if instrument is em.Inst.TRECS:
                        data_label = header.get('DATALAB')
                        if data_label is not None and '-CAL' in data_label:
                            result = ObservationIntentType.CALIBRATION
                    elif instrument is em.Inst.CIRPASS:
                        data_label = header.get('DATALAB')
                        if data_label is not None and '-CAL' in data_label:
                            result = ObservationIntentType.CALIBRATION
                        else:
                            obs_type = header.get('OBS_TYPE')
                            if obs_type in cal_values:
                                result = ObservationIntentType.CALIBRATION
                            else:
                                result = ObservationIntentType.SCIENCE
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
    if instrument is em.Inst.PHOENIX:
        result = _get_phoenix_obs_type(header)
    elif instrument is em.Inst.HOKUPAA:
        # DB - 01-18-19 Use IMAGETYP as the keyword
        result = header.get('IMAGETYP')
    elif instrument is em.Inst.CIRPASS:
        # DB - 06-03-19
        temp = header.get('OBJECT')
        # if ‘dome’ or ‘flat’ or ‘twilight’ or ‘gcal’ in
        #       object string:  obstype = FLAT
        # if ‘argon’ or ‘arc’ in string: obstype = ARC
        # if ‘dark’ in string: obstype = DARK
        # Otherwise assume obstype = OBJECT
        result = 'OBJECT'
        for ii in ['dome', 'flat', 'twilight', 'gcal']:
            if ii in temp:
                result = 'FLAT'
                break
        for ii in ['argon', 'arc']:
            if ii in temp:
                result = 'ARC'
                break
        if 'dark' in temp:
            result = 'DARK'
    return result


def get_proposal_id(header):
    """
    Determine the Proposal ID.

    :param header:  The FITS header for the current extension.
    :return: The proposal id from Gemini JSON metadata, or None if not found.
    """
    return em.om.get('program_id')


def get_ra(header):
    """
    Get the right ascension. Rely on the JSON metadata, because it's all in
    the same units (degrees).

    :param header:  The FITS header for the current extension.
    :return: ra, or None if not found.
    """
    instrument = _get_instrument()
    if instrument is em.Inst.HOKUPAA:
        ra, dec = _get_sky_coord(header, 'RA', 'DEC')
        result = ra
    elif instrument is em.Inst.OSCIR:
        ra, dec = _get_sky_coord(header, 'RA_TEL', 'DEC_TEL')
        result = ra
    elif instrument in [em.Inst.BHROS, em.Inst.NIFS, em.Inst.TEXES]:
        # bHROS, TEXES: ra/dec not in json
        # DB - 05-03-19 - NIFS: ra/dec not reliable in json
        result = header.get('RA')
    elif instrument is em.Inst.CIRPASS:
        # DB - 06-03-19 - Must use FITS header info for most WCS info
        ra, dec = _get_sky_coord(header, 'TEL_RA', 'TEL_DEC')
        result = ra
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


def get_target_type(uri):
    """
    Calculate the Target TargetType

    :param header:  The FITS header for the current extension.
    :return: The Target TargetType, or None if not found.
    """
    em.om.reset_index(uri)
    result = TargetType.FIELD
    spectroscopy = em.om.get('spectroscopy')
    instrument = _get_instrument()
    if spectroscopy or instrument is em.Inst.TEXES:
        result = TargetType.OBJECT
    return result


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

        # DB - 06-03-19 json ut_date_time Gemini is based on DATE keyword.
        # Better to use MJD header keyword value directly.

        time_val = header.get('MJD')
    else:
        # DB - 07-03-19 - TEXES json ut_date_time is correct.
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


def _get_data_label():
    return em.om.get('data_label')


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
    data_label = _get_data_label()
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


def _is_gmos_mask(header):
    result = False
    instrument = _get_instrument()
    if instrument in [em.Inst.GMOSS, em.Inst.GMOSN, em.Inst.GMOS]:
        obs_type = _get_obs_type(header)
        if obs_type == 'MASK':
            result = True
    return result


def accumulate_fits_bp(bp, file_id, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug('Begin accumulate_fits_bp for {}.'.format(file_id))
    em.get_obs_metadata(file_id)
    bp.set('Observation.intent', 'get_obs_intent(header)')
    bp.set('Observation.type', 'get_obs_type(header)')
    bp.set('Observation.metaRelease', 'get_meta_release(header)')
    bp.set('Observation.target.type', 'get_target_type(uri)')
    bp.set('Observation.target.moving', 'get_target_moving(header)')
    bp.set('Observation.proposal.id', 'get_proposal_id(header)')

    telescope = em.om.get('telescope')
    if telescope is not None:
        if 'North' in telescope:
            x, y, z = ac.get_location(19.823806, -155.46906, 4213.0)
        else:
            x, y, z = ac.get_location(-30.240750, -70.736693, 2722.0)
        bp.set('Observation.telescope.geoLocationX', x)
        bp.set('Observation.telescope.geoLocationY', y)
        bp.set('Observation.telescope.geoLocationZ', z)

    bp.set('Plane.productID', file_id)
    bp.set('Plane.dataProductType', 'get_data_product_type(header)')
    bp.set('Plane.calibrationLevel', 'get_calibration_level(uri)')
    bp.set('Plane.metaRelease', 'get_meta_release(header)')
    bp.set('Plane.dataRelease', 'get_data_release(header)')

    bp.set('Plane.provenance.name', 'Gemini Observatory Data')
    bp.set('Plane.provenance.project', 'Gemini Archive')
    # Add IMAGESWV for GRACES
    bp.add_fits_attribute('Plane.provenance.producer', 'IMAGESWV')
    bp.set_default('Plane.provenance.producer', 'Gemini Observatory')
    instrument = _get_instrument()
    if instrument is not em.Inst.TEXES:
        data_label = _get_data_label()
        bp.set('Plane.provenance.reference',
               'http://archive.gemini.edu/searchform/{}'.format(data_label))

    bp.set('Artifact.productType', 'get_art_product_type(header)')
    bp.set('Artifact.contentChecksum', 'md5:{}'.format(em.om.get('data_md5')))
    bp.set('Artifact.contentLength', em.om.get('data_size'))
    bp.set('Artifact.contentType', 'application/fits')
    # always see the metadata, see the data only when it's public
    bp.set('Artifact.releaseType', 'data')
    bp.set('Artifact.uri', uri)

    if instrument is em.Inst.CIRPASS:
        bp.set_default('Observation.telescope.name', 'Gemini-South')
    mode = em.om.get('mode')
    if not (instrument in [em.Inst.GPI, em.Inst.PHOENIX, em.Inst.HOKUPAA,
                           em.Inst.OSCIR, em.Inst.BHROS] or
            (instrument is em.Inst.GRACES and mode is not None and
             mode != 'imaging')):
        bp.configure_position_axes((1, 2))

    if instrument is em.Inst.FLAMINGOS:
        # DB 27-05-19
        # Flamingos, you actually want to use the EQUINOX value, not the
        # EPOCH.   And I think EQUINOX header value is usually 2000.0, even
        # for the example GS-CAL20020620-15-0462 02jun20.0462 with
        # RA_TEL = “UNAVAILABLE”.  For Gemini the assumption is that the
        # RA/Dec values in the headers are always based on the position of
        # the equinox given at the time specified by the EQUINOX keyword value.
        bp.clear('Chunk.position.equinox')
        bp.add_fits_attribute('Chunk.position.equinox', 'EQUINOX')
        
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
    if 'product_id' in kwargs:
        current_product_id = kwargs['product_id']

    # processed files
    if is_composite(headers):
        logging.info('{} is a Composite Observation.'.format(
            observation.observation_id))
        observation = _update_composite(observation)

    instrument = em.Inst(observation.instrument.name)

    if instrument in [em.Inst.MICHELLE, em.Inst.GNIRS]:
        # DB 16-04-19
        # The more important issue with this and other files is that they
        # contain no image extensions.  The file is downloadable from
        # the Gemini archive but their only content is the primary
        # header.   There is no pixel data.  Test for the existence of a
        # FITS extension and skip processing of a michelle file if
        # there isn’t one

        # DB 18-04-19
        #
        # For the last GNIRS file (NAXIS=0)  skip the file if it doesn’t have
        # an extension.

        if len(headers) == 1:
            logging.warning(
                '{}: no image data for {}. Cannot build an '
                'observation.'.format(instrument, observation.observation_id))
            return None

    try:
        for p in observation.planes:
            if current_product_id != observation.planes[p].product_id:
                continue

            for a in observation.planes[p].artifacts:

                caom_name = ec.CaomName(observation.planes[p].artifacts[a].uri)
                em.om.reset_index(caom_name.uri)
                processed = em.gofr.is_processed(caom_name.file_name)
                if (instrument in
                        [em.Inst.MICHELLE, em.Inst.TRECS, em.Inst.GNIRS]):
                    # Michelle is a retired visitor instrument.
                    # Spatial WCS info is in primary header. There
                    # are a variable number of FITS extensions
                    # defined by primary keyword NUMEXT; assume the
                    # same spatial WCS for each for now - it differs
                    # only slightly because of telescope 'chopping'
                    # and 'nodding' used in acquisition. DB - 01-18-19
                    #
                    # DB - 01-18-19 - GNIRS has no WCS info in extension; use
                    # primary header
                    _update_position_from_zeroth_header(
                        observation.planes[p].artifacts[a], headers,
                        instrument, observation.observation_id)

                for part in observation.planes[p].artifacts[a].parts:

                    if part == '2' and instrument is em.Inst.GPI:
                        # GPI data sets have two extensions. First is science
                        # image (with WCS), second is data quality for each
                        # pixel (no WCS).
                        logging.info(
                            'GPI: Setting chunks to None for part {} for {}'.format(
                                part, observation.observation_id))
                        observation.planes[p].artifacts[a].parts[part].chunks \
                            = TypedList(Chunk, )
                        continue
                    for c in observation.planes[p].artifacts[a].parts[
                        part].chunks:
                        header = headers[int(part)]

                        # energy WCS
                        filter_name = get_filter_name(
                            headers[0], header, observation.observation_id,
                            instrument)
                        if _reset_energy(observation.type, p, instrument,
                                         filter_name):
                            c.energy = None
                            c.energy_axis = None
                        else:
                            if instrument is em.Inst.NIRI:
                                _update_chunk_energy_niri(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument is em.Inst.GPI:
                                _update_chunk_energy_gpi(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument is em.Inst.F2:
                                _update_chunk_energy_f2(
                                    c, headers[0],
                                    observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument is em.Inst.GSAOI:
                                _update_chunk_energy_general(
                                    c, instrument, [DataProductType.IMAGE],
                                    observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument is em.Inst.HRWFS:
                                _update_chunk_energy_hrwfs(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument is em.Inst.NICI:
                                _update_chunk_energy_nici(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument is em.Inst.TRECS:
                                _update_chunk_energy_trecs(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument is em.Inst.MICHELLE:
                                _update_chunk_energy_michelle(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument is em.Inst.NIFS:
                                _update_chunk_energy_nifs(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument is em.Inst.GNIRS:
                                _update_chunk_energy_gnirs(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                            elif instrument is em.Inst.PHOENIX:
                                _update_chunk_energy_phoenix(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument is em.Inst.FLAMINGOS:
                                _update_chunk_energy_flamingos(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id,
                                    filter_name)
                                ignore, obs_type = _get_flamingos_mode(header)
                                if (obs_type is not None and
                                        observation.type is None):
                                    observation.type = obs_type
                            elif instrument is em.Inst.HOKUPAA:
                                _update_chunk_energy_hokupaa(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument is em.Inst.OSCIR:
                                _update_chunk_energy_oscir(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name)
                            elif instrument is em.Inst.BHROS:
                                _update_chunk_energy_bhros(
                                    c, header,
                                    observation.planes[p].data_product_type,
                                    observation.observation_id)
                            elif instrument is em.Inst.GRACES:
                                _update_chunk_energy_graces(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id)
                            elif instrument in [em.Inst.GMOS, em.Inst.GMOSN,
                                                em.Inst.GMOSS]:
                                _update_chunk_energy_gmos(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id, filter_name,
                                    instrument)
                            elif instrument is em.Inst.CIRPASS:
                                _update_chunk_energy_cirpass(
                                    c, observation.planes[p].data_product_type,
                                    observation.observation_id)
                            elif instrument is em.Inst.TEXES:
                                _update_chunk_energy_texes(
                                    c, headers[0],
                                    observation.planes[p].data_product_type,
                                    observation.observation_id)

                        # position WCS
                        mode = em.om.get('mode')
                        if _reset_position(headers, instrument):
                            logging.debug(
                                'Setting Spatial WCS to None for {}'.format(
                                    observation.observation_id))
                            c.position_axis_2 = None
                            c.position_axis_1 = None
                            c.position = None
                        else:
                            if (instrument in [em.Inst.PHOENIX, em.Inst.HOKUPAA,
                                               em.Inst.OSCIR] or
                                    (instrument is em.Inst.GRACES and
                                     mode is not None and mode != 'imaging')):
                                _update_chunk_position(
                                    c, header, instrument, int(part),
                                    observation.observation_id)
                            elif instrument in [em.Inst.BHROS, em.Inst.CIRPASS,
                                                em.Inst.TEXES]:
                                _update_chunk_position(
                                    c, headers[0], instrument,
                                    int(part), observation.observation_id)
                            elif instrument is em.Inst.NIFS:
                                # DB - 01-18-19 - NIFS has no WCS info in
                                # extension; use primary header
                                #
                                # DB - 04-03-19 - NIFS spatial WCS info in the
                                # header has way too large a FOV so hard-code
                                # this to the instrument's tiny 3" x 3" FOV.

                                n_axis1 = headers[-1]['NAXIS1']
                                n_axis2 = headers[-1]['NAXIS2']
                                _update_chunk_position(
                                    c, headers[0], instrument,
                                    int(part), observation.observation_id,
                                    n_axis1, n_axis2)
                            elif instrument is em.Inst.GPI:
                                _update_chunk_position(
                                    c, headers[1], instrument, int(part),
                                    observation.observation_id)
                                if part == '1':
                                    # equinox information only available from
                                    # 0th header
                                    c.position.equinox = headers[0].get('TRKEQUIN')
                            elif instrument is em.Inst.FLAMINGOS:
                                _update_chunk_position_flamingos(
                                    c, header, observation.observation_id)

                        # time WCS
                        if instrument is em.Inst.F2:
                            _update_chunk_time_f2(c, observation.observation_id)

                if isinstance(observation, CompositeObservation):
                    cc.update_plane_provenance(observation.planes[p],
                                               headers[1:],
                                               'IMCMB', COLLECTION,
                                               _repair_provenance_value,
                                               observation.observation_id)

                if ((processed or
                     isinstance(observation, CompositeObservation) or
                     instrument is em.Inst.TEXES)
                        and 'jpg' not in caom_name.file_name):
                    # not the preview artifact
                    if observation.planes[p].provenance is not None:
                        observation.planes[p].provenance.reference = \
                            'http://archive.gemini.edu/searchform/' \
                            'filepre={}'.format(caom_name.file_name)

            program = em.get_pi_metadata(observation.proposal.id)
            if program is not None:
                observation.proposal.pi_name = program['pi_name']
                observation.proposal.title = program['title']

        if isinstance(observation, CompositeObservation):
            cc.update_observation_members(observation)

    except Exception as e:
        logging.error('Error {} for {} instrument {}'.format(
            e, observation.observation_id, instrument))
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

# DB 21-05-19
# For lack of info on the slit for these odd observations use the
# same resolution for the 4-pixel slit:  i.e. add “‘f6-cam’: 610.0”
# for the J grism lookup (of which there may only be two cases)
# and “‘f6-cam’: 780.0" for the K grism.  I don’t see any similar
# types of observations for other grisms that are not ‘blank’
# observations and so don’t have energy WCS.

NIRI_RESOLVING_POWER = {
    'J': {
        'f6-2pix': 770.0,
        'f6-4pix': 610.0,
        'f6-cam': 610.0,
        'f6-6pix': 460.0,
        'f6-2pixBl': 770.0,
        'f6-4pixBl': 650.0,
        'f6-6pixBl': 480.0,
        'f32-4pix': 1000.0,
        'f32-6pix': 620.0,  # f32-7pix
        'f32-9pix': 450.0  # f32-10pix
    },
    'H': {
        'f6-2pix': 1650.0,
        'f6-4pix': 825.0,
        'f6-6pix': 520.0,
        'f6-2pixBl': 1650.0,
        'f6-4pixBl': 940.0,
        'f6-6pixBl': 550.0,
        'f32-4pix': 880.0,
        'f32-6pix': 630.0,  # f32-7pix
        'f32-9pix': 500.0  # f32-10pix
    },
    'L': {
        'f6-2pix': 1100.0,
        'f6-4pix': 690.0,
        'f6-6pix': 460.0,
        'f6-2pixBl': 1100.0,
        'f6-4pixBl': 770.0,
        'f6-6pixBl': 490.0,
    },
    'M': {
        'f6-2pix': 1100.0,
        'f6-4pix': 770.0,
        'f6-6pix': 460.0
    },
    'K': {
        'f6-2pix': 1300.0,
        'f6-4pix': 780.0,
        'f6-cam': 780.0,
        'f6-6pix': 520.0,
        'f6-2pixBl': 1300.0,
        'f6-4pixBl': 780.0,
        'f6-6pixBl': 520.0,
        'f32-4pix': 1280.0,
        'f32-6pix': 775.0,  # f32-7pix
        'f32-9pix': 570.0  # f32-10pix
    }
}

# NIRI.Hgrismf32-G5228


def _update_chunk_energy_niri(chunk, data_product_type, obs_id, filter_name):
    """NIRI-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_niri')
    mc.check_param(chunk, Chunk)

    reset_energy = False

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

    # DB - 08-09-19 - NIRI data with ‘pinhole’:  ignore energy WCS for
    # these for any spectroscopic observations.  They are used for image
    # quality tests and instrument alignment.  For imaging observations
    # these shouldn’t impact anything (I hope).  e.g. does this one work?
    # GN-2005A-DD-12-12-017

    # DB - 15-04-19
    # For r_ratio:
    # https://www.gemini.edu/sciops/instruments/niri/spectroscopy/grisms,
    # parenthetically states “w. either slits” for the f/6 K section.  So
    # duplicate lines 1268-1270 of main_app.py with “Bl” appended after
    # ‘pix’.

    if 'Jcon(112)_G0235' in filter_name:
        # DB - 01-04-19 The G0235 filter is listed as ‘damaged’ on the Gemini
        # NIRI filters web site:
        # https://www.gemini.edu/sciops/instruments/niri/imaging/filters.
        # Not enough info is given there for SVO to add this filter to their
        # system.  Hardcode a central wavelength of 1.1232 microns and a
        # FWHM of 0.0092 microns
        filter_md = FilterMetadata('NIRI')
        filter_md.central_wl = 1.1232
        filter_md.bandpass = 0.0092
    elif 'Msort' in filter_name or 'Mgrism' in filter_name:
        # DB - 15-04-19
        # The NIRI blocking filter page,
        # https://www.gemini.edu/sciops/instruments/niri/spectroscopy/blocking-
        # filters, gives no transmission data for this filter for the SVO
        # to use.  Will need to hard-code it to central wavelength of
        # (4.4+6)/2 and width of 6-4.4 microns.

        filter_md = FilterMetadata('NIRI')
        filter_md.set_central_wl(6., 4.4)
        filter_md.set_bandpass(6., 4.4)
    else:
        filter_md = em.get_filter_metadata(em.Inst.NIRI, filter_name)
        if filter_md is None:
            raise mc.CadcException('{}: mystery filter {} for {}'.format(
                em.Inst.NIRI, filter_name, obs_id))

    filter_name = em.om.get('filter_name')
    if data_product_type == DataProductType.IMAGE:
        logging.debug('NIRI: SpectralWCS imaging for {}.'.format(obs_id))
        fm = filter_md
        fm.adjust_resolving_power()
    elif data_product_type == DataProductType.SPECTRUM:
        logging.debug('NIRI: SpectralWCS spectroscopy for {}.'.format(
            obs_id))
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
            if (bandpass_name in NIRI_RESOLVING_POWER and
                    f_ratio in NIRI_RESOLVING_POWER[bandpass_name]):
                fm.resolving_power = \
                    NIRI_RESOLVING_POWER[bandpass_name][f_ratio]
            elif 'pinhole' in f_ratio:
                logging.info(
                    'Pinhole. Setting energy to None for {}'.format(obs_id))
                reset_energy = True
            else:
                raise mc.CadcException(
                    'NIRI: Mystery bandpass name {} or f_ratio {} for '
                    '{}.'.format(bandpass_name, f_ratio, obs_id))
        else:
            raise mc.CadcException(
                'NIRI: Mystery disperser value {} for {}'.format(disperser,
                                                                 obs_id))
    else:
        raise mc.CadcException(
            'NIRI: Do not understand mode {} for {}'.format(
                data_product_type, obs_id))

    if reset_energy:
        chunk.energy_axis = None
        chunk.energy = None
    else:
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


def _update_chunk_energy_f2(chunk, header, data_product_type, obs_id,
                            filter_name):
    """F2-specific chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_f2')
    mc.check_param(chunk, Chunk)

    # DB - 02-05-19
    # For F2 use SVO filter service for bandpass or ‘delta’ for images and
    # spectroscopy.  Treat images as for other instruments but…
    #
    # if gemini_md[‘spectroscopy’] == ‘true’:
    #  ref_wl = gemini_md[‘central_wavelength’]   or  GRWLEN header value
    #  grism = gemini_md[‘disperser’]  or  GRISM header value
    #  if gemini_md[‘mode’] == ‘LS’:  # long-slit
    #    slit_width = MASKNAME header value, e.g. ‘4pix-slit’, but need only ‘4’
    #    use the table I sent you with slit/grism values to determine
    #    average resolution R
    #
    #  elif gemini_md[‘mode’] == ‘MOS’:  # Multi-object
    #    slit_width = 2   # no way to determine slit widths used in
    #    custom mask, so assume 2
    #    use the table I sent you with slit/grism values to determine
    #    average resolution R
    #
    #  else:
    #    fail because there shouldn’t be any other spectroscopy mode

    # DB 09-04-19 - Ignore energy when the grism is in the header but
    # object value of “COVER CLOSED” so is another type of calibration
    # exposure apparently.

    reset_energy = False
    object_value = em.om.get('object')
    if 'COVER CLOSED' in object_value or 'Undefined' in filter_name:
        # DB 30-04-19
        # Flamingos ‘Undefined’ filter:  no spectral WCS
        reset_energy = True
    else:
        filter_md = em.get_filter_metadata(em.Inst.F2, filter_name)
        if data_product_type == DataProductType.IMAGE:
            logging.debug('SpectralWCS: F2 imaging mode for {}.'.format(obs_id))
            fm = filter_md
        elif data_product_type == DataProductType.SPECTRUM:
            logging.debug(
                'SpectralWCS: F2 LS|Spectroscopy mode for {}.'.format(obs_id))
            fp_mask = header.get('MASKNAME')
            mode = em.om.get('mode')
            slit_width = None
            if mode == 'LS':
                slit_width = fp_mask[0]
            fm = FilterMetadata()
            fm.central_wl = filter_md.central_wl
            fm.bandpass = filter_md.bandpass
            grism_name = header.get('GRISM')
            logging.debug(
                'F2: grism name is {} fp_mask is {} for {}'.format(grism_name,
                                                                   fp_mask,
                                                                   obs_id))
            # lookup values from
            # https://www.gemini.edu/sciops/instruments/flamingos2/spectroscopy/longslit-spectroscopy
            lookup = {'1': [1300.0, 3600.0],
                      '2': [900.0, 2800.0],
                      '3': [600.0, 1600.0],
                      '4': [350.0, 1300.0],
                      '6': [130.0, 1000.0],
                      '8': [100.0, 750.0]}
            if slit_width is None or slit_width not in lookup:
                # DB 02-04-19
                # For F2 at line 1409 of main_app.py set slit_width = ‘2’ as
                # a default of slit_width[0] is not a numeric value
                slit_width = '2'
            if grism_name.startswith('R3K_'):
                fm.resolving_power = lookup[slit_width][1]
            else:
                fm.resolving_power = lookup[slit_width][0]
        else:
            raise mc.CadcException(
                'F2: Do not understand DataProductType {} for {}'.format(
                    data_product_type, obs_id))

    if reset_energy:
        logging.info(
            'Setting spectral WCs to none for {} instrument {}'.format(
                obs_id, em.Inst.F2))
        chunk.energy_axis = None
        chunk.energy = None
    else:
        _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_f2')


def _update_chunk_energy_general(chunk, instrument, allowable_types,
                                 data_product_type, obs_id, filter_name):
    """General chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_{}'.format(instrument))
    mc.check_param(chunk, Chunk)

    if data_product_type in allowable_types:
        logging.debug('{} Spectral WCS {} mode for {}.'.format(instrument,
                                                               data_product_type,
                                                               obs_id))
        filter_md = em.get_filter_metadata(instrument, filter_name)
        if filter_md is None:
            raise mc.CadcException(
                '{}: mystery filter {}'.format(instrument, filter_name))
        _build_chunk_energy(chunk, filter_name, filter_md)
    else:
        raise mc.CadcException(
            '{} no Spectral WCS support when DataProductType {} for {}'.format(
                instrument, data_product_type, obs_id))

    logging.debug('End _update_chunk_energy_{}'.format(instrument))


def _update_chunk_energy_hrwfs(chunk, data_product_type, obs_id, filter_name):
    """General chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_hrwfs')
    mc.check_param(chunk, Chunk)

    # DB 27-05-19
    # e.g. GS-CAL20020322-7-0003 2002mar22_0055, filter1=‘neutral’,
    # filter2=‘open’, so treat similarly to “open” for GMOS since it’s a
    # similar CCD:  maybe 0.35 - 1.0 microns.

    if data_product_type == DataProductType.IMAGE:
        logging.debug(
            'hrwfs Spectral WCS {} mode for {}.'.format(data_product_type,
                                                        obs_id))
        if 'open' in filter_name or 'neutral' in filter_name:
            w_min = 0.35
            w_max = 1.0
            filter_md = FilterMetadata()
            filter_md.set_bandpass(w_max, w_min)
            filter_md.set_central_wl(w_max, w_min)
            filter_md.set_resolving_power(w_max, w_min)
            # DB 27-05-19
            # bandpassName likely best to set to NULL
            filter_name = ''
        else:
            filter_md = em.get_filter_metadata(em.Inst.HRWFS, filter_name)
            temp = []
            for ii in filter_name.split('+'):
                temp.append(ii.split('_')[0])
            filter_name = '+'.join(i for i in temp)
        _build_chunk_energy(chunk, filter_name, filter_md)
    else:
        raise mc.CadcException(
            'hrwfs no Spectral WCS support when DataProductType {} for '
            '{}'.format(data_product_type, obs_id))

    logging.debug('End _update_chunk_energy_hrwfs')


# 0 - central
# 1 - lower
# 2 - upper
#
# units are microns
#
# select filter_id,wavelength_central,
# wavelength_lower,wavelength_upper from gsa..gsa_filters where
# instrument="michelle";
michelle = {
    'F103B10': [10.282990, 9.875740, 10.690240],
    'F105B53': [10.500000, 7.850000, 13.150000],
    'F112B21': [11.299610, 10.236910, 12.362310],
    'F116B9': [11.687970, 11.254020, 12.121920],
    'F125B9': [12.493580, 12.097480, 12.889680],
    'F128B2': [12.800000, 12.650000, 12.950000],
    'F14SA': [10.000000, 1.000000, 14.000000],
    'F161L': [20.500000, 16.100000, 25.000000],
    'F185B9B': [18.113860, 17.385410, 18.842310],
    'F198B27': [20.696080, 17.945430, 23.446730],
    'F209B42L': [20.900000, 16.500000, 25.300000],
    'F209B42S': [20.900000, 16.500000, 25.300000],
    'F66LA': [10.000000, 6.600000, 99.900000],
    'F66LB': [10.000000, 6.600000, 99.900000],
    'F79B10': [7.715890, 7.392540, 8.039240],
    'F88B10': [8.821500, 8.469400, 9.173600],
    'F97B10': [9.688450, 9.253750, 10.123150],
    'QBlock': [(16.1 + 25) / 2, 16.1, 25.],
    'NBlock': [(14 + 6.8) / 2, 6.8, 14.],
}


def _update_chunk_energy_michelle(chunk, data_product_type, obs_id,
                                  filter_name):
    """General chunk-level Energy WCS construction."""
    logging.debug('Begin _update_chunk_energy_michelle')
    mc.check_param(chunk, Chunk)

    # DB - 01-03-19
    #
    # no resolution for imaging mode
    #
    # Michelle spectroscopy:
    # json ‘disperser’, e.g. LowN, MedN1, MedN2, LowQ, Echelle
    # json ‘focal_plane_mask’, e.g. 2_pixels (need the ‘2’)
    #
    # Disp.    R    final R as function of slit width
    # --------------------------------------------------
    # LowN    200    R x 2/slit width
    # LowQ    110    R x 3/slit width
    # MedN1    1000    R x 2/slit width
    # MedN2    3000    R x 2/slit width
    # Echelle    30000    R x 2/slit width (very approximate)

    # DB - 04-04-19
    # The only solution for the michelle datasets (another ‘visitor’
    # instrument) would be to use filter info in the gsa..gsa_filters
    # table (as you do for PHOENIX).  No info to pass on to SVO folks to
    # add more filters.   Code would have to handle this case of two
    # filters with overlapping bandpasses.

    # DB - 15-04-19
    # Michelle:  From a hidden page of Michelle filters,
    # http://www.gemini.edu/sciops/instruments/michelle/imaging/filters,
    # hard-code QBlock and NBlock filters using values in that table
    # (ignoring the greater-than symbols)?  i.e. NBlock has central
    # bandpass of  (14+6.8)/2 microns and bandpass of 14-6.8 microns.
    # Ditto for QBlock:  (16.1+25)/2 and 25-16.1.  Grid_T should be
    # ignored for bandpass calculations but it would be good to keep it
    # in the filter name so in this case it would be F125B9 + Grid_T

    # 0 = R
    # 1 = ratio for slit width
    lookup = {'LowN': [200, 2.0],
              'LowQ': [110, 3.0],
              'MedN1': [1000, 2.0],
              'MedN2': [3000, 2.0],
              'Echelle': [30000, 2.0]
              }

    filter_md = em.get_filter_metadata(em.Inst.MICHELLE, filter_name)
    if filter_md is None:  # means filter_name not found
        # DB 09-04-19 - Use 100 microns for the initial max for michelle.
        w_max, w_min = _multiple_filter_lookup(filter_name, michelle,
                                               obs_id, em.Inst.MICHELLE,
                                               wl_max=100)
        filter_md = FilterMetadata()
        filter_md.set_bandpass(w_max, w_min)
        filter_md.set_central_wl(w_max, w_min)
    if data_product_type == DataProductType.SPECTRUM:
        logging.debug(
            'michelle: Spectral WCS spectrum for {}.'.format(obs_id))
        if data_product_type == DataProductType.SPECTRUM:
            disperser = em.om.get('disperser')
            focal_plane_mask = em.om.get('focal_plane_mask')
            slit_width = float(focal_plane_mask.split('_')[0])
            if disperser not in lookup:
                raise mc.CadcException(
                    'michelle: Mystery disperser {} for {}'.format(disperser,
                                                                   obs_id))
            filter_md.resolving_power = \
                lookup[disperser][0] * lookup[disperser][1] / slit_width
    elif data_product_type == DataProductType.IMAGE:
        logging.debug(
            'michelle: Spectral WCS imaging mode for {}.'.format(obs_id))
    else:
        raise mc.CadcException(
            'michelle: no Spectral WCS support when DataProductType {} for '
            '{}'.format(data_product_type, obs_id))

    # use the json value for bandpass_name value - it's representative of
    # multiple filters
    filter_name = em.om.get('filter_name')
    _build_chunk_energy(chunk, filter_name, filter_md)
    logging.debug('End _update_chunk_energy_michelle')


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
    if filter_name == 'Block':
        # DB 04-04-19
        # If one of the NICI filters is ‘Block’ then energy WCS should be
        # ignored for that extension.
        chunk.energy = None
        chunk.energy_axis = None
    else:
        filter_md = em.get_filter_metadata(em.Inst.NICI, filter_name)

        if data_product_type == DataProductType.IMAGE:
            logging.debug(
                'NICI: SpectralWCS imaging mode for {}.'.format(obs_id))
            if filter_md is None:  # means filter_name not found
                w_max, w_min = _multiple_filter_lookup(filter_name, NICI,
                                                       obs_id, em.Inst.NICI)
                fm = FilterMetadata()
                fm.set_bandpass(w_max, w_min)
                fm.set_central_wl(w_max, w_min)
                fm.set_resolving_power(w_max, w_min)
            else:
                fm = filter_md

            temp = em.om.get('filter_name')
            # NICI has two different bandpass names (most of the time) in
            # its two chunks.  Pat says in this case nothing will be put in
            # the bandpass name for the plane.  Add code to combine the two
            # chunk bandpass names to create a plane bandpass name only for
            # this instrument
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
    if filter_md is None:
        raise mc.CadcException(
            '{}: Mystery filter {}'.format(em.Inst.TRECS, filter_name))
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
            'TReCS: Do not understand mode {} for {}'.format(data_product_type,
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

    # DB - 04-04-19
    # NIFS gratings can be tuned (rotated) for slightly different
    # wavelengths.  But these lines should be added to the NIFS lookup:
    #
    # ‘J’: {‘JH’: [1.25, 1.15, 1.33, 6040.0, 49.6]}
    # ‘H’: {‘HK’: [1.65, 1.49, 1.80, 5290.0, 56.8]}
    #
    # i.e. the same values as the other J/H grating entries except for
    # the filter names.  What’s important are the 2nd and 3rd numbers
    # that are used to determine the bandpass since you use the json
    # ‘central_wavelength’ to establish that value.  That central
    # wavelength combined with the bandpass should set different
    # upper/lower wavelength limits for NIFS observations with the same
    # grating/filter combination.
    #
    # For the NIFS JH filter observations with disperser K those appear
    # to be observations that they use to measure noise or dark signal or
    # perhaps point to a very bright star for an acquisition observation.
    # There are other combinations I see, e.g. ZJ with K.  Those should
    # skip energy WCS calculation.  i.e. any filter/disperser combination
    # not in the lookup.   Also, I see a ‘Blocked’ filter for NIFS.
    # Those should also skip energy WCS.

    nifs_lookup = {'Z': {'ZJ': [1.05, 0.94, 1.15, 4990.0, 60.1]},
                   'J': {'ZJ': [1.25, 1.15, 1.33, 6040.0, 49.6],
                         'JH': [1.25, 1.15, 1.33, 6040.0, 49.6]},
                   'H': {'JH': [1.65, 1.49, 1.80, 5290.0, 56.8],
                         'HK': [1.65, 1.49, 1.80, 5290.0, 56.8]},
                   'K': {'HK': [2.20, 1.99, 2.40, 5290.0, 56.7]},
                   'K_Short': {'HK': [2.20, 1.98, 2.41, 5290.0, 56.7]},
                   'K_Long': {'HK': [2.20, 1.98, 2.41, 5290.0, 56.7]}}

    if 'INVALID' in filter_name or 'Blocked' in filter_name:
        fm = None
    else:
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
                    logging.info(
                        'NIFS: No energy. filter_name {} with disperser {}'
                        ' for {}'.format(filter_name, grating, obs_id))
                    fm = None
            else:
                logging.info(
                    'NIFS: No energy. grating {} for {}'.format(
                        grating, obs_id))
                fm = None
        elif data_product_type == DataProductType.IMAGE:
            logging.debug('NIFS: imaging for {}.'.format(obs_id))
            # DB - 01-03-19
            # NIFS images should just use the standard imaging procedure for
            # resolution (central_wavelength/bandpass).
        else:
            raise mc.CadcException(
                'NIFS: DataProductType {} for {}'.format(
                    data_product_type, obs_id))

    if fm is None:
        chunk.energy_axis = None
        chunk.energy = None
    else:
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

    # long slit mode information source:
    # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy
    # 0 - min wavelength
    # 1 - max wavelength
    # 2 - 'short' camera resolution
    # 3 - 'long' camera resolution
    # 4 - Since November 2012 and for the cross-dispersed mode with
    # the 2 pix wide slit only resolving powers are somewhat lower, as
    # follows: X-1400; J-1400, H-1400; K-1300, for 'short' camera
    # resolution

    # cross-dispersion information:
    # Change xd_mode to include grating ID (e.g. ‘32’ or ‘10’ or ‘110’)
    # and two more configurations:
    # xd_mode = {‘SB+SXD+32’: [0.9, 2.5, 1800.0],
    # ‘LB+LXD+10’: [0.9, 2.5, 1800.0],
    # ‘LB+SXD+10’: [1.2, 2.5, 1800.0],
    # ‘SB+SXD+110’:[0.9, 2.5, 5400.0],
    # ‘LB+SXD+110’:[0.9, 2.5, 5400.0],
    # ‘SB+LXD+32’:{0.9, 2.5, 5400.0]}

    # DB - 05-03-19
    # Disregard central wavelength provided by Gemini.  e.g. if a long-slit
    # observation with the K filter in the beam has a central_wavelength
    # setting of (for example) 2.3 microns then using the bounding
    # wavelengths in the arrays to define the wavelength coverage will
    # result in the upper wavelength being outside the filter bandpass.
    # Instead use the average of the lower/upper wavelength ranges
    # for each configuration in long_slit_mode and xd_mode as the central
    # wavelength.  e.g. for K long-slit observations the central wavelength
    # would be (1.91+2.49)/2.0 or 2.2 microns.  The way it is now with a
    # central wavelength of 2.3 microns we calculate wavelength limits of
    # 2.01 to 2.59 microns but wavelengths beyond 2.49 microns don’t make it
    # past the filter.

    # DB - 03-04-19

    # GNIRS observation GN-2017A-Q-44-25-031 has a json ‘mode’ value of
    # ‘imaging’ despite having a disperser value of ‘32_mm&SXD’ indicating
    # that it’s a spectrum (as you can also see from the preview). Means
    # forgetting about relying on json ‘mode’ but looking to see if the
    # json disperser value = ‘MIRROR’.

    # DB 08-09-19 - observations like GN-2010B-SV-142-655-012 with null
    # 'disperser'.  Some kind of rare calibration exposure - there are
    # only on the order of 100 in the archive - or just an error in the
    # header.  The spectra look ‘dispersed’ to me.  Just ignore energy
    # WCS for these.

    # DB 24-04-19
    # ND = neutral density and so any ND* filter can be ignored as it
    # shouldn’t affect transmission band.  Likely observing a bright
    # target and they need to reduce light by a factor of 100X.

    gnirs_lookup = {'10': {'X': [1.03, 1.17, 570, 2100],
                           'J': [1.17, 1.37, 570, 1600],
                           'H': [1.47, 1.80, 570, 1700],
                           'K': [1.91, 2.49, 570, 1700],
                           'L': [2.80, 4.20, 570, 1800],
                           'M': [4.40, 6.00, 570, 1200],
                           'LB+LXD': [0.9, 2.5, 1800],
                           'LB+SXD': [1.2, 2.5, 1800]},
                    '32': {'X': [1.03, 1.17, 1700, 5100, 1400],
                           'J': [1.17, 1.37, 1600, 4800, 1400],
                           'H': [1.49, 1.80, 1700, 5100, 1400],
                           'K': [1.91, 2.49, 1700, 5100, 1300],
                           'L': [2.80, 4.20, 1800, 5400, 1800],
                           'M': [4.40, 6.00, 1240, 3700, 1240],
                           'SB+SXD': [0.9, 2.5, 1800, 1800],
                           'SB+LXD': [0.9, 2.5, 5400, 5400],
                           'LB+LXD': [0.9, 2.5, 5400, 5400],
                           'LB+SXD': [1.2, 2.5, 1800, 1800]},
                    '111': {'X': [1.03, 1.17, 6600, 17800],
                            'J': [1.17, 1.37, 7200, 17000],
                            'H': [1.49, 1.80, 5900, 17800],
                            'K': [1.91, 2.49, 5900, 17800],
                            'L': [2.80, 4.20, 6400, 19000],
                            'M': [4.40, 6.00, 4300, 12800],
                            'SB+SXD': [0.9, 2.5, 5400],
                            'LB+SXD': [0.9, 2.5, 5400],
                            'LB+LXD': [0.9, 2.5, 17000]}}

    reset_energy = False
    if 'Dark' in filter_name:
        # 'Dark' test obs is GN-2013B-Q-93-147-036
        # DB 13-05-19
        # GNIRS “Dark” trumps the “L” filter, so no energy.
        reset_energy = True
        logging.warning(
            'GNIRS: filter is {}. No spectral WCS for {}.'.format(
                filter_name, obs_id))
    else:
        fm = FilterMetadata('GNIRS')
        if data_product_type == DataProductType.SPECTRUM:
            logging.debug(
                'gnirs: SpectralWCS Spectroscopy mode for {}.'.format(obs_id))
            disperser = em.om.get('disperser')
            if disperser is None:
                logging.info('No disperser. No energy for {}'.format(obs_id))
                reset_energy = True
            else:
                grating = disperser.split('_')[0]
                # 'UNKNOWN' in grating test obs GS-CAL20040924-6-006
                # 'ENG -' in grating test obs GN-CAL20130813-22-010
                if ('UNKNOWN' in grating or 'Moving' in grating or
                        'ENG' in grating):
                    # DB 23-04-19 - GNIRS grating UNKNOWN:  no energy
                    # DB 18-04-19
                    # If grating is moving then the observation is almost
                    # certainly garbage so ignore energy WCS.
                    # DB 13-05-19
                    # No energy for the invalid “ENG - 170000" grating entry.
                    logging.warning('grating is {}. No energy for {}'.format(
                        grating, obs_id))
                    reset_energy = True
                else:
                    if grating not in gnirs_lookup:
                        raise mc.CadcException(
                            'GNIRS: Mystery grating {} for {}'.format(grating, obs_id))

                    camera = em.om.get('camera')
                    focal_plane_mask = em.om.get('focal_plane_mask')
                    if focal_plane_mask is None or 'arcsec' not in focal_plane_mask:
                        # DB 24-04-19
                        # Assume slit width of 1.0 for GNIRS observation
                        # without a focal plane mask value.
                        slit_width = 1.0
                    else:
                        slit_width = float(focal_plane_mask.split('arcsec')[0])

                    if 'XD' in disperser:
                        logging.debug(
                            'gnirs: cross dispersed mode for {}.'.format(obs_id))
                        # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/
                        # crossdispersed-xd-spectroscopy/xd-prisms
                        # 0 = lower
                        # 1 = upper
                        # 2 = spectral resolution with 2-pix wide slit
                        # DB - 04-03-19 - Change the last number in each row to
                        # 1800.0 since the resolving power is the same for all 3
                        # cases
                        #
                        # Add line to find grating ID (from long-slit code):
                        # grating = disperser.split(‘_’)[0]
                        #
                        # Then change ‘lookup’ to include grating.
                        #
                        # I can’t find any other combinations (e.g. ‘LB+LXD+32’)
                        # but no guarantee that I won’t have to add another line
                        # or two if we see failures.   Wavelength coverage isn’t
                        # correct for the R=5400 entries because only about
                        # 1/3rd of the full band pass is covered but in bits
                        # and pieces that we can’t identify in raw image.

                        lookup = None
                        coverage = disperser[-3:]
                        if camera.startswith('Short'):
                            lookup = '{}+{}'.format('SB', coverage)
                        elif camera.startswith('Long'):
                            lookup = '{}+{}'.format('LB', coverage)

                        if camera.startswith('Long'):
                            slit_table_value = 0.1
                            lookup_index = 2
                        elif camera.startswith('Short'):
                            slit_table_value = 0.3
                            lookup_index = 2
                        else:
                            raise mc.CadcException(
                                'GNIRS: Mystery camera definition {} for {}'.format(
                                    camera,
                                                                                    obs_id))
                    else:
                        logging.debug('gnirs: long slit mode for {}.'.format(
                            obs_id))
                        bandpass = filter_name[0]
                        lookup = bandpass
                        if camera.startswith('Long'):
                            slit_table_value = 0.1
                            lookup_index = 3
                        elif camera.startswith('Short'):
                            date_time = ac.get_datetime(em.om.get('ut_datetime'))
                            if date_time > ac.get_datetime('2012-11-01T00:00:00'):
                                slit_table_value = 0.3
                                if grating == '32':
                                    lookup_index = 4
                                else:
                                    lookup_index = 3
                            else:
                                slit_table_value = 0.3
                                lookup_index = 2
                        else:
                            raise mc.CadcException(
                                'GNIRS: Mystery camera definition {} for {}'.format(
                                    camera, obs_id))

                    if lookup not in gnirs_lookup[grating]:
                        raise mc.CadcException(
                            'GNIRS: Mystery lookup {} for grating {}, obs {}'.format(
                                lookup, grating, obs_id))
                    bounds = gnirs_lookup[grating][lookup]
                    fm.set_bandpass(bounds[1], bounds[0])
                    fm.resolving_power = slit_table_value * bounds[
                        lookup_index] / slit_width
                    fm.set_central_wl(bounds[1], bounds[0])
        elif data_product_type == DataProductType.IMAGE:
            logging.debug('gnirs: SpectralWCS imaging mode for {}.'.format(obs_id))
            # https://www.gemini.edu/sciops/instruments/gnirs/imaging
            # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/orderblocking-acq-nd-filters

            # DB 23-04-19
            # GNIRS filter “L”: add lower/upper bandpass info for the L and M
            # filter on this page,
            # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/orderblocking-acq-nd-filters
            # to ‘imaging’ dictionary.  (In ‘imaging’ I think the K order
            # blocking filter lower bandpass should be 1.91).
            # Despite Gemini’s footnote stating that L and M filters are NOT
            # used for acquisition images.

            # 0 - min wavelength
            # 1 - max wavelength
            imaging = {'Y': [0.97, 1.07],
                       'J': [1.17, 1.33],
                       'J order blocking)': [1.17, 1.37],
                       'H': [1.49, 1.80],
                       'K': [2.03, 2.37],
                       'K order blocking': [1.91, 2.49],
                       'H2': [2.105, 2.137],
                       'PAH': [3.27, 3.32],
                       'X': [1.03, 1.17],
                       'XD': [0.9, 2.56],
                       'M': [4.4, 6.0],
                       'L': [2.8, 4.2]}

            # DB 30-04-19
            # GN-CAL20180607-27-001:  ignore ‘PupilViewer’ if possible.
            filter_name = filter_name.replace('PupilViewer', '')
            filter_name = filter_name.strip('+')

            if (len(filter_name) == 1 or filter_name == 'H2' or
                    filter_name == 'PAH' or filter_name == 'XD'):
                bandpass = filter_name
            else:
                bandpass = filter_name[0]

            if bandpass in imaging:
                bounds = imaging[bandpass]
            else:
                raise mc.CadcException(
                    'GNIRS: Unexpected filter_name {} for {}'.format(
                        filter_name, obs_id))

            fm.set_central_wl(bounds[1], bounds[0])
            fm.set_bandpass(bounds[1], bounds[0])
        else:
            raise mc.CadcException(
                'GNIRS: Unexpected DataProductType {} for {}'.format(
                    data_product_type, obs_id))

    if reset_energy:
        chunk.energy_axis = None
        chunk.energy = None
    else:
        _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_gnirs')


# select filter_id, wavelength_central, wavelength_lower, wavelength_upper
# from gsa..gsa_filters where instrument = 'PHOENIX'
# 0 - central
# 1 - lower
# 2 - upper
#
# units are microns
#
# dict with the filter wheels stripped from the names as returned by query
#
# DB 30-04-19
# A new filter not previously available.  Add info from here,
# https://www.noao.edu/kpno/phoenix/filters.html

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
           'L2870': [3.490500, 3.436000, 3.545000],
           '9440': [1.058500, 1.053000, 1.064000]}


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


def _update_chunk_energy_flamingos(chunk, data_product_type, obs_id,
                                   filter_name):
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
        if fm is None:
            raise mc.CadcException(
                'Flamingos: Mystery filter {} for {}'.format(filter_name,
                                                             obs_id))

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
                      'Kp': [2.12, 0.41],
                      'H+K notch': [1.8, 0.7],
                      'methane low': [1.56, 0.0120],
                      'methane high': [1.71, 0.0120],
                      'FeII': [1.65, 0.0170],
                      'HeI': [2.06, 0.0030],
                      '1-0 S(1) H2': [2.12, 0.0023],
                      'H Br(gamma)': [2.166, 0.0150],
                      'K-continuum': [2.26, 0.0060],
                      'CO': [2.29, 0.0020]
                      }

    # DB 27-05-19
    # The bandpasses for the 7 bottom entries have to be corrected because of
    # the misleading column heading.  The bandpasses for these 7 are in
    # Angstroms, not microns.  So divide the values by 10,000 to convert to
    # microns. 170 Å = 0.017 microns for the FeII filter.
    # H2/23 = 1-0 S(1) H2 - bandpass should be 0.0023
    # HKnotch = H+K notch
    # 1.56/120 = methane low  - bandpass should be 0.0120
    # 1.71/120 = methane high - ditto
    # LowFlx is likely the J+CO combo in line 14 of the filter table in the
    # web page:  effectively a shutter so no energy in this case.
    # And I think FeI/17 is likely meant to be the FeII filter (just going by
    # the “17” in the bandpass.  I’m not aware of an IR Fe I filter…

    filter_name_repair = {'H2/23': '1-0 S(1) H2',
                          'HKnotch': 'H+K notch',
                          '1.56/120': 'methane low',
                          '1.71/120': 'methane high',
                          'FeI/17': 'FeII'}
    reset_energy = False

    if 'LowFlx' in filter_name:
        reset_energy = True
    else:
        if filter_name in filter_name_repair:
            filter_name = filter_name_repair[filter_name]
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

    if reset_energy:
        chunk.energy_axis = None
        chunk.energy = None
    else:
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

    # DB 23-04-19
    # Split on whitespace and name the filters S_8.8 and IHW_(17-19)
    # respectively.  And change the keys of the lookup table to these
    # values as well.Randomly looking at observations these are other
    # filter values:  S_12.5 (-11774), S_11.7 (-22125), S_7.9 (-63440),
    # and N_wide (-1154).   The latter instead of the lookup “N”.  So I’m
    # guessing the first six key values in oscir_lookup should be
    # S_7.9, S_8.8, S_9.8, S_11.7 and S_12.5.   Then N_wide, IHW_(17-19).
    #
    # 0 - central wavelenth
    # 1 - bandpass
    # units are microns
    oscir_lookup = {'S_7.9': [7.91, 0.755],
                    'S_8.8': [8.81, 0.871],
                    'S_9.8': [9.80, 0.952],
                    'S_10.3': [10.27, 1.103],
                    'S_11.7': [11.70, 1.110],
                    'S_12.5': [12.49, 1.156],
                    'N_wide': [10.75, 5.230],
                    'IHW_(17-19)': [18.17, 1.651],
                    'Q3': [20.8, 1.650]}

    temp = filter_name
    filter_name = temp.split()[0]
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

    # DB - 04-04-19

    # GG455&g filter is a case where both the GG455 order-sorting filter
    # and the g filter are combined (which is likely quite common).  The
    # former is in the ‘lookup’ on line 2343.  The g filter is g_G0301
    # which is in the SVO.  My original ‘svo’ script handled the
    # order-sorting filters as part of that subroutine.  I don’t think
    # that works with the current _update_chunk_energy_gmos.  The o-s
    # filter will change the effective wavelength coverage of the g
    # filter.

    # this means some filter metadata will come from SVO, and some is
    # hard-coded here

    # DB 18-04-19
    #
    # Lots of missing header data for the GMOS observation
    # GN-2005B-Q-60-11-011 so ignore energy WCS

    reset_energy = False

    filter_md = None
    if 'open' not in filter_name and 'No_Value' not in filter_name:
        filter_md = em.get_filter_metadata(instrument, filter_name)
    w_max = 10.0
    w_min = 0.0
    for ii in filter_name.split('+'):
        if 'Hartmann' in ii:
            continue
        elif ii in lookup:
            wl_max = lookup[ii][1]
            wl_min = lookup[ii][0]
            if wl_max < w_max:
                w_max = wl_max
            if wl_min > w_min:
                w_min = wl_min

    if filter_md is not None:
        # mingle the SVO lookup values with the hard-coded lookup values
        # from here
        svo_min = (2 * filter_md.central_wl - filter_md.bandpass) / 2
        svo_max = filter_md.bandpass + svo_min
        if svo_max < w_max:
            w_max = svo_max
        if svo_min > w_min:
            w_min = svo_min

    filter_md = FilterMetadata()
    filter_md.set_bandpass(w_max, w_min)
    filter_md.set_central_wl(w_max, w_min)
    filter_md.set_resolving_power(w_max, w_min)

    if data_product_type == DataProductType.SPECTRUM:
        logging.debug('{}: SpectralWCS spectroscopy for {}.'.format(
            instrument, obs_id))
        if math.isclose(filter_md.central_wl, 0.0):
            logging.info(
                '{}: no spectral wcs, central wavelength is {} for {}'.format(
                    instrument, filter_md.central_wl, obs_id))
            return
        disperser = em.om.get('disperser')
        # 'unknown' in disperser test obs is GN-2004B-Q-30-15-002
        # 'OLDMIRROR' in disperser test obs is GN-CAL20020329-2-025
        # 'No Value' in disperser test obs is GN-2013B-SV-152-120-001
        if (disperser is None or 'unknown' in disperser or
                'OLDMIRROR' in disperser or 'No Value' in disperser):
            # DB 23-04-19
            # GMOS observation with OLDMIRROR looks to be a direct image
            # with the pinhole mask in the beam.  But Gemini json info is
            # incorrect since it claims it’s a spectrum.  No energy.
            #
            # GMOS observation GN-CAL20020501-3-000 with unknown
            # disperser.  No energy.  Just the observation number 000
            # flags this one as unusual.
            #
            # DB 13-05-19
            # “No Value” for filter?  Yes, no energy WCS.  Lots of other bogus
            # values in header as well.
            logging.warning('{}: disperser is {}, no energy.'.format(
                instrument, disperser))
            reset_energy = True
        else:
            fm = FilterMetadata()
            fm.central_wl = filter_md.central_wl
            fm.bandpass = filter_md.bandpass
            if disperser == 'B12000':
                # DB 16-04-19
                # B12000 must be a Gemini typo since observation
                # GN-2006B-Q-39-100-003 has B1200 for the disperser.
                disperser = 'B1200'
            elif disperser in ['B600-', 'B600+-G5323', '\'B600']:
                # DB 23-04-19
                # GMOS observation GS-CAL20030130-1-002 with B600- grating.
                # The ‘-’ must be a typo.  The observation preceding this
                # one, GS-CAL20030130-1-001, has the B600 grating.
                # DB 13-05-19
                # GS-2013A-Q-91-194-003 S20130517S0092
                # Yes, disperser value should be B600.  Typo using +- instead
                # of the usual underscore when someone entered the info in a
                # config file perhaps?
                disperser = 'B600'
            if disperser in GMOS_RESOLVING_POWER:
                fm.resolving_power = GMOS_RESOLVING_POWER[disperser]
            else:
                raise mc.CadcException(
                    'gmos: mystery disperser {} for {}'.format(disperser, obs_id))
    elif data_product_type == DataProductType.IMAGE:
        logging.debug('{}: SpectralWCS imaging for {}.'.format(
            instrument, obs_id))
        fm = filter_md
    else:
        raise mc.CadcException(
            '{}: mystery data product type {} for {}'.format(
                instrument, data_product_type, obs_id))
    if reset_energy:
        chunk.energy_axis = None
        chunk.energy = None
    else:
        _build_chunk_energy(chunk, filter_name, fm)
    logging.debug('End _update_chunk_energy_gmos')


def _update_chunk_energy_cirpass(chunk, data_product_type, obs_id):
    # DB - 06-03-19
    # Energy WCS:
    #
    # Can’t do anything better than fixing lower/upper wavelength bounds to
    # 1.0 and 1.67 microns and resolving power of 3200.  NAXIS=1 as we’ve
    # done for other spectra like this.
    #
    # Data type should always be spectrum despite json giving ‘imaging’.
    logging.debug('Begin _update_chunk_energy_cirpass')
    if data_product_type == DataProductType.SPECTRUM:
        logging.debug(
            'cirpass: SpectralWCS spectral mode for {}.'.format(obs_id))
        fm = FilterMetadata()
        fm.set_central_wl(1.0, 1.67)
        fm.set_bandpass(1.0, 1.67)
        fm.resolving_power = 3200.0
    else:
        raise mc.CadcException(
            'cirpass: mystery data product type {} for {}'.format(
                data_product_type, obs_id))
    _build_chunk_energy(chunk, '', fm)
    logging.debug('End _update_chunk_energy_cirpass')


def _update_chunk_energy_texes(chunk, header, data_product_type, obs_id):
    # DB - 07-03-19
    # TEXES Spectroscopy
    #
    # Some special code will be needed for datalabels/planes.  There are no
    # datalabels in the FITS header.  json metadata (limited) must be
    # obtained with URL like
    # https://archive.gemini.edu/jsonsummary/canonical/filepre=TX20170321_flt.2507.fits.
    # Use TX20170321_flt.2507 as datalabel.  But NOTE:  *raw.2507.fits and
    # *red.2507.fits are two planes of the same observation. I’d suggest we
    # use ‘*raw*’ as the datalabel and ‘*red*’ or ‘*raw*’ as the appropriate
    # product ID’s for the science observations.  The ‘flt’ observations do
    # not have a ‘red’ plane.  The json document contains ‘filename’ if
    # that’s helpful at all.  The ‘red’ files do not exist for all ‘raw’
    # files.
    #
    #
    # Header OBSTYPE appears to be correct; not json obs_type.
    #
    # No previews are generated by Gemini
    #
    # Energy WCS:
    #
    # Central wavelength is given by 10,000/header(WAVENO0).  I have to do
    # some more investigation to see if we can determine wavelength coverage
    # (i.e. see if I can identify the echelle/echelon info rom the header - I
    # don’t think so).  For now use 0.25 microns as the fixed FWHM bandpass.
    logging.debug('Begin _update_chunk_energy_texes')
    if data_product_type == DataProductType.SPECTRUM:
        logging.debug(
            'texes: SpectralWCS spectral mode for {}.'.format(obs_id))
        fm = FilterMetadata('TEXES')
        fm.central_wl = 10000/header.get('WAVENO0')
        fm.bandpass = 0.25
    else:
        # data_type/observing mode is always spectroscopy
        raise mc.CadcException(
            'texes: mystery data product type {} for {}'.format(
                data_product_type, obs_id))
    _build_chunk_energy(chunk, '', fm)
    logging.debug('End _update_chunk_energy_texes')


def _reset_energy(observation_type, data_label, instrument, filter_name):
    """
    Return True if there should be no energy WCS information created at
    the chunk level.

    :param observation_type from the parent Observation instance.
    :param data_label str for useful logging information only.
    :param filter_name other information useful for knowing whether to
        reset energy
    """
    result = False
    om_filter_name = em.om.get('filter_name')

    if ((observation_type is not None and
         ((observation_type == 'DARK') or
          (instrument in [em.Inst.GMOS, em.Inst.GMOSN, em.Inst.GMOSS] and
           observation_type in ['BIAS', 'MASK']))) or
            (om_filter_name is not None and ('blank' in om_filter_name or
                                             'Blank' in om_filter_name)) or
            (filter_name is not None and ('unknown' in filter_name or
                                          filter_name == ''))):
        logging.info(
            'No chunk energy for {} obs type {} filter name {}'.format(
                data_label, observation_type, om_filter_name))

        # 'unknown' in filter_name test obs is GN-2004B-Q-30-15-002
        # DB 23-04-19 - GN-2004B-Q-30-15-002: no energy
        # GMOS GS-2005A-Q-26-12-001.  Lots of missing metadata, including
        # release date so no energy (filter_name == '')
        result = True
    elif instrument is em.Inst.GNIRS and 'Moving' in filter_name:
        # DB 16-04-19
        # Energy WCS should be ignored for ‘Moving’ since we don’t know
        # what might have been in the light path when the exposure was
        # actually being taken.
        result = True
    elif instrument is em.Inst.TRECS and filter_name == '':
        # DB 23-04-19
        # GMOS GS-2005A-Q-26-12-001.  Lots of missing metadata, including
        # release date so no energy.  Ditto for TReCS GS-CAL20041206-6-007.
        result = True
    elif (instrument is em.Inst.MICHELLE and
          ('Blank' in filter_name or 'No Value' in filter_name)):
        # 'No Value' in filter_name test obs GN-2005A-C-14-45-002
        # DB 09-04-19 - “Blank-B” -> no energy.  It is a bias exposure
        # apparently and those shouldn't have energy WCS.
        # DB 23-04-19 - “No Value” -> no energy
        result = True
    elif (instrument is em.Inst.NIRI and
          ('INVALID' in om_filter_name or filter_name == '')):
        # filter_name == '', test obs is GN-CAL20050301-17-001
        # DB 18-04-19
        #
        # Gemini archive shows WaveBand=INVALID - ignore energy WCS
        # DB 23-04-19 - ‘no value’ -> no energy

        # DB 23-04-19 - Looks like NIRI observation GN-CAL20020623-1-011 is
        # skipping the ‘invalid’ values in FILTER1/2 headers.  So energy
        # should likely be skipped.  In this case the json value is
        # INVALID&INVALID.
        result = True
    elif (instrument is em.Inst.GSAOI and (
            'Unknown' in filter_name or 'Blocked' in filter_name)):
        # DB 24-04-19
        # ‘Unknown+Blocked2’ filter, no spectral WCS.
        result = True
    return result


def _reset_position(headers, instrument):
    """
    Return True if there should be no spatial WCS information created at
    the chunk level.
    """
    result = False
    types = em.om.get('types')
    ra = get_ra(headers[0])
    if 'AZEL_TARGET' in types and ra is None:
        # DB - 02-04-19 - Az-El coordinate frame likely means the telescope
        # was parked or at least not tracking so spatial information is
        # irrelevant.

        # DB - 09-04-19 - AZEL_TARGET should likely be checked for all
        # datasets, and means the spatial WCS should be ignored. since this
        # generally means the telescope is not tracking and so spatial WCS
        # info isn’t relevant since the position is changing with time.
        result = True
    elif instrument in [em.Inst.NIFS, em.Inst.PHOENIX]:
        # DB - 08-04-19 - json ra/dec values are null for
        # the file with things set to -9999.  Ignore
        # spatial WCS for these cases.

        # get the values from JSON directly, because the
        # function uses header values, which are set to
        # unlikely defaults

        # DB 30-04-19
        # Looks like many relatively recent PHOENIX files have no RA/Dec
        # values in the header and so will have no spatial WCS.
        # Base this decision on json null values.  But looking at all
        # of the PHOENIX data from 2016 until 3 December 2017 it looks
        # like json ra/dec values are either null or 0.0 for all.
        # In both cases spatial WCS should be ignored.  (It will be
        # very difficult for users to find anything of interest in
        # these datasets other than searching by free-form target
        # names…)  PHOENIX returned as a visitor instrument in May 2016
        # after about 5 years away.

        ra = em.om.get('ra')
        dec = em.om.get('dec')
        if ra is None and dec is None:
            result = True
        elif (ra is not None and math.isclose(ra, 0.0) and
              dec is not None and  math.isclose(dec, 0.0)):
            result = True
    elif _is_gmos_mask(headers[0]):
        # DB - 04-03-19
        # Another type of GMOS-N/S dataset to archive.
        # Mask images.   json observation_type = “MASK”.
        # These have no WCS info at all, although I guess
        # json ut_date_time could be used as the start date
        # with null exposure time. These would have only
        # instrument, obstype, datatype (spectrum) and
        # product type (AUXILIARY) set.
        result = True
    elif instrument is em.Inst.GRACES:
        # DB 23-04-19
        # Ignore spatial WCS for the GRACES dataset with EPOCH=0.0.  Not
        # important for a bias. For GMOS we skip spatial WCS for biases
        # (and maybe for some other instruments).

        # DB 24-04-19
        # GRACES:  you can ignore spatial WCS for flats if RA/Dec are not
        # available.   Ditto for GNIRS darks.

        # DB 30-04-19
        # Ignore spatial WCS for any GRACES arcs without RA/Dec values.

        ra = em.om.get('ra')
        dec = em.om.get('dec')
        if ra is None and dec is None:
            obstype = get_obs_type(headers[0])
            if obstype in ['BIAS', 'FLAT', 'ARC']:
                result = True
    elif instrument is em.Inst.FLAMINGOS:
        ra_tel = headers[0].get('RA_TEL')
        if ra_tel == 'Unavailable':
            result = True
    return result


def _multiple_filter_lookup(filter_name, lookup, obs_id, instrument, wl_max=None):
    w_max = 10.0 if wl_max is None else wl_max
    w_min = 0.0
    for ii in filter_name.split('+'):
        if ii in lookup:
            wl_max = lookup[ii][2]
            wl_min = lookup[ii][1]
        else:
            msg = '{}: Unprepared for filter {} from {}'.format(
                instrument, ii, obs_id)
            if (instrument == em.Inst.MICHELLE and
                    (ii.startswith('I') or (ii == 'Grid_T'))):
                logging.info(msg)
                continue
            else:
                raise mc.CadcException(msg)
        if wl_max < w_max:
            w_max = wl_max
        if wl_min > w_min:
            w_min = wl_min
    return w_max, w_min


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
    if instrument in [em.Inst.GRACES, em.Inst.TEXES]:
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
    # MICHELLE
    # DB - 08-04-19
    # Use header FILTERA and FILTERB values to determine the filter
    # bandpass
    #
    if instrument not in [em.Inst.NIRI, em.Inst.NICI, em.Inst.MICHELLE]:
        filter_name = em.om.get('filter_name')
        #
        # DB 24-04-19
        # ND = neutral density and so any ND* filter can be ignored as it
        # shouldn’t affect transmission band.  Likely observing a bright
        # Other instruments occasionally have ND filters in the beam.
        if filter_name is not None:
            filter_name = filter_name.replace('&', '+')
            temp = filter_name.split('+')
            for fn in temp:
                if fn.startswith('ND'):
                    filter_name = filter_name.replace(fn, '')
            filter_name = filter_name.strip('+')
    if (filter_name is None or
            (filter_name is not None and len(filter_name.strip()) == 0)):
        # DB - 04-02-19 - strip out anything with 'pupil' as it doesn't affect
        # energy transmission
        #
        # 08-04-19 - MICHELLE
        # Note that FILTERB sometimes has a value “Clear_B” and that should
        # be ignored.  Not sure if “Clear_A” is possible but maybe best to
        # allow for it.
        #
        # DB 24-04-19
        # Other instruments occasionally have ND filters in the beam.
        filters2ignore = ['open', 'invalid', 'pupil', 'clear', 'nd']
        lookups = ['FILTER']
        if instrument is em.Inst.PHOENIX:
            lookups = ['CVF_POS', 'FILT_POS']
            search_header = header
        elif instrument is em.Inst.NICI:
            search_header = header
        else:
            search_header = primary_header
        for key in search_header.keys():
            for lookup in lookups:
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


def _update_position_from_zeroth_header(artifact, headers, instrument, obs_id):
    """Make the 0th header spatial WCS the WCS for all the
    chunks."""
    primary_header = headers[0]

    # naxis values are only available from extensions

    primary_header['NAXIS1'] = headers[-1].get('NAXIS1')
    primary_header['NAXIS2'] = headers[-1].get('NAXIS2')

    primary_chunk = Chunk()
    types = em.om.get('types')
    ra = get_ra(headers[0])
    ra_json = em.om.get('ra')
    dec_json = em.om.get('dec')
    if (('AZEL_TARGET' in types and ra is None) or
            (instrument is em.Inst.GNIRS and ra_json is None and
             dec_json is None)):
        # DB - 02-04-19 - Az-El coordinate frame likely means the telescope
        # was parked or at least not tracking so spatial information is
        # irrelevant.

        # DB - 09-04-19 - AZEL_TARGET should likely be checked for all
        # datasets, and means the spatial WCS should be ignored. since this
        # generally means the telescope is not tracking and so spatial WCS
        # info isn’t relevant since the position is changing with time.

        # DB 24-04-19
        # Ignore spatial WCS if RA/Dec are not available for GNIRS.

        logging.info(
            '{}: Spatial WCS is None for {}'.format(instrument, obs_id))
    else:
        wcs_parser = WcsParser(primary_header, obs_id, 0)
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


def _update_chunk_position(chunk, header, instrument, extension, obs_id,
                           n_axis1=None, n_axis2=None):
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

    # TEXES - DB - 07-03-19
    # Use header RA and DEC for crval1/2.  Use a fixed 5" x 5" FOV.
    # So NAXIS=1 and cd11=cd22 = 5.0/3600.0

    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CRVAL1'] = get_ra(header)
    header['CRVAL2'] = get_dec(header)
    if instrument not in [em.Inst.GPI, em.Inst.CIRPASS]:
        header['CDELT1'] = RADIUS_LOOKUP[instrument]
        header['CDELT2'] = RADIUS_LOOKUP[instrument]
        header['CROTA1'] = 0.0
    if instrument not in [em.Inst.OSCIR, em.Inst.GPI]:
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
    if instrument is em.Inst.NIFS:
        # DB 05-03-19 - persist NAXIS values for NIFS
        header['NAXIS1'] = n_axis1
        header['NAXIS2'] = n_axis2
    if instrument is em.Inst.CIRPASS:
        # So perhaps try:
        #     NAXIS1 = 33
        #     NAXIS2 = 15
        header['NAXIS1'] = 33
        header['NAXIS2'] = 15
        # LENS_SCL determines the scale/lenslet:  0.36 or 0.25 (arcseconds
        # per lens)
        #     if 0.36 then FOV is 13.0" x 4.7" (RA and Dec)
        #     if 0.25 then FOV is 9.3" x 3.5"
        lens_scl = header.get('LENS_SCL')
        if lens_scl == '0.36':
            header['CDELT1'] = 13.0 / 3600.0
            header['CDELT2'] = 4.7 / 3600.0
        else:
            header['CDELT1'] = 9.3 / 3600.0
            header['CDELT2'] = 3.5 / 3600.0
    header['CRPIX1'] = get_crpix1(header)
    header['CRPIX2'] = get_crpix2(header)
    header['CD1_1'] = get_cd11(header)
    header['CD1_2'] = 0.0
    header['CD2_1'] = 0.0
    header['CD2_2'] = get_cd22(header)

    if instrument is em.Inst.BHROS:
        header['EQUINOX'] = float(header.get('TRKEQUIN'))

    if instrument is em.Inst.PHOENIX:
        temp = header.get('EQUINOX')
        if temp is None or math.isclose(temp, 0.0):
            header['EQUINOX'] = header.get('EPOCH')

    wcs_parser = WcsParser(header, obs_id, extension)
    if chunk is None:
        chunk = Chunk()
    wcs_parser.augment_position(chunk)
    chunk.position_axis_1 = 1
    chunk.position_axis_2 = 2

    if instrument is em.Inst.OSCIR:
        chunk.position.coordsys = header.get('FRAMEPA')
        chunk.position.equinox = float(header.get('EQUINOX'))
    elif instrument is em.Inst.BHROS:
        chunk.position.coordsys = header.get('TRKFRAME')
    elif instrument is em.Inst.GPI:
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


def _update_chunk_time_gmos(chunk, obs_id):
    """"""
    logging.debug('Begin _update_chunk_time_gmos {}'.format(obs_id))
    mc.check_param(chunk, Chunk)
    if chunk.time is not None and chunk.time.axis is not None:
        mjd_start = get_time_function_val(None)
        start = RefCoord(0.5, mjd_start.value)
        end = RefCoord(1.5, mjd_start.value)
        chunk.time.exposure = 0.0
        chunk.time.axis.range = CoordRange1D(start, end)
        chunk.time.axis.function = None
    logging.debug('End _update_chunk_time_gmos {}'.format(obs_id))


def _update_composite(obs):
    comp_obs = change_to_composite(obs)
    return comp_obs


def _repair_provenance_value(imcmb_value, obs_id):
    """There are several naming patterns in the provenance for
    processed files. Try to extract meaningful raw file names.
    Return None if the encountered pattern is unexpected."""
    if 'N' in imcmb_value:
        temp = 'N' + imcmb_value.split('N', 1)[1]
    elif 'S' in imcmb_value:
        temp = 'S' + imcmb_value.split('S', 1)[1]
    elif '$' in imcmb_value:
        temp = imcmb_value.split('$', 1)[1]
    else:
        logging.warning(
            'Unrecognized IMCMB value {}'.format(imcmb_value))
        return None, None

    if '_' in temp:
        temp1 = temp.split('_')[0]
    elif '.fits' in temp:
        temp1 = temp.split('.fits')[0]
    elif '[SCI' in temp:
        temp1 = temp.split('[SCI')[0]
    else:
        logging.warning(
            'Failure to repair {} for {}'.format(temp, obs_id))
        return None, None

    prov_file_id = temp1[:14]
    prov_obs_id = em.gofr.get_obs_id(prov_file_id)
    if prov_obs_id is None:
        return None, prov_file_id
    return prov_obs_id[0], prov_file_id


def _build_blueprints(uris):
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
            file_id = GemName.remove_extensions(ec.CaomName(uri).file_name)
            accumulate_fits_bp(blueprint, file_id, uri)
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


def is_composite(headers):
    """All the logic to determine if a file name is part of a
    CompositeObservation, in one marvelous function."""
    result = False

    # look in the last header - IMCMB keywords are not in the zero'th header
    header = headers[-1]
    for ii in header:
        if ii.startswith('IMCMB'):
            result = True
            break
    return result


def change_to_composite(observation):
    """For the case where a SimpleObservation needs to become a
    CompositeObservation."""
    return CompositeObservation(observation.collection,
                                observation.observation_id,
                                Algorithm('composite'),
                                observation.sequence_number,
                                observation.intent,
                                observation.type,
                                observation.proposal,
                                observation.telescope,
                                observation.instrument,
                                observation.target,
                                observation.meta_release,
                                observation.planes,
                                observation.environment,
                                observation.target_position)


def main_app2():
    args = get_gen_proc_arg_parser().parse_args()
    if em.gofr is None:
        em.gofr = GemObsFileRelationship('/app/data/from_paul.txt')
    try:
        uris = _get_uris(args)
        blueprints = _build_blueprints(uris)
        gen_proc(args, blueprints)
    except Exception as e:
        logging.error('Failed {} execution for {}.'.format(APPLICATION, args))
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)

    logging.debug('Done {} processing.'.format(APPLICATION))
