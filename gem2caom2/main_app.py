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

# TODO - none of this is really true anymore, so figure out what needs
to be captured here for documentation ..... ;)

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
import re
import sys
import traceback

from astropy import units
from astropy.coordinates import SkyCoord

from caom2 import Observation, ObservationIntentType, DataProductType
from caom2 import CalibrationLevel, TargetType, ProductType, Chunk, Axis
from caom2 import SpectralWCS, CoordAxis1D, RefCoord, Instrument
from caom2 import TypedList, CoordRange1D, DerivedObservation
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2utils import WcsParser
from caom2pipe import manage_composable as mc
from caom2pipe import caom_composable as cc
from caom2pipe import astro_composable as ac

import gem2caom2.external_metadata as em
import gem2caom2.obs_file_relationship as ofr
from gem2caom2.gem_name import GemName, COLLECTION
from gem2caom2.svofps import FilterMetadata
from gem2caom2.builder import get_instrument, GemObsIDBuilder
from gem2caom2 import instruments


__all__ = ['gem_main_app', 'to_caom2', 'update', 'APPLICATION']

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
RADIUS_LOOKUP = {
    em.Inst.GPI: 2.8 / 3600.0,  # units are arcseconds
    em.Inst.GRACES: 1.2 / 3600.0,
    em.Inst.PHOENIX: 5.0 / 3600.0,
    em.Inst.OSCIR: 0.0890 / 3600.0,
    em.Inst.HOKUPAA: 4.0 / 3600.0,
    em.Inst.BHROS: 0.9 / 3600.0,
    em.Inst.NIFS: 3.0 / 3600.0,
    em.Inst.TEXES: 5.0 / 3600.0,
}


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
    instrument = get_instrument()
    if (
        (
            reduction is not None
            and (('PROCESSED' in reduction) or ('PREPARED' in reduction))
        )
        or (
            instrument is em.Inst.TEXES
            and ('_red' in uri.lower() or '_sum' in uri.lower())
        )
        or (
            instrument is em.Inst.PHOENIX
            and mc.CaomName(uri.lower()).file_id.startswith('p')
        )
        or (
            instrument is em.Inst.OSCIR
            and mc.CaomName(uri.lower()).file_id.startswith('r')
        )
    ):
        # DB 23-02-21
        # The best thing to do with OSCIR 'r' files is to add them as a
        # second cal level 2 plane and use the same metadata as the
        # unprocessed plane.
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
        f'obs type is {obs_type} obs class is {obs_class} for '
        f'{_get_data_label()}'
    )
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
            instrument = get_instrument()
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
    instrument = get_instrument()
    if instrument is em.Inst.HOKUPAA:
        result = _get_pix_scale(header)
    elif instrument in [em.Inst.OSCIR, em.Inst.TEXES]:
        result = RADIUS_LOOKUP[instrument]
    elif instrument in [em.Inst.GPI, em.Inst.NIFS]:
        # DB - 05-03-19 - NIFS needs a division by NAXIS1/2 for the
        # cdelta1/2 calculations.
        result = RADIUS_LOOKUP[instrument] / header.get('NAXIS1')
    elif instrument is em.Inst.CIRPASS:
        # DB - 06-03-19
        # FOV is fixed at two possible values and has no bearing on NAXIS1/2
        # values. See
        # http://www.gemini.edu/sciops/instruments/cirpass/cirpassIFU.html.
        # 499 lenslets cover the FOV: about 33 along one axis and 15 along the
        # other.
        #
        # LENS_SCL determines the scale/lenslet:  0.36 or 0.25
        # (arcseconds per lens)
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
    instrument = get_instrument()
    if instrument is em.Inst.HOKUPAA:
        result = _get_pix_scale(header)
    elif instrument in [em.Inst.OSCIR, em.Inst.TEXES]:
        result = RADIUS_LOOKUP[instrument]
    elif instrument in [em.Inst.GPI, em.Inst.NIFS]:
        # DB - 05-03-19 - NIFS needs a division by NAXIS1/2 for the
        # cdelta1/2 calculations.
        result = RADIUS_LOOKUP[instrument] / header.get('NAXIS2')
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
    instrument = get_instrument()
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
    instrument = get_instrument()
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
    instrument = get_instrument()
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
    elif instrument in [em.Inst.ALOPEKE, em.Inst.ZORRO]:
        # DB 31-08-20
        # Both cameras are used for speckle imaging.  Datasets consist of
        # cubes of 256 x 256 x 1000 images.  i.e. 1000 short exposures << 1
        # second long with 256 x 256 pixel images (or smaller images if the CCD
        # is binned).
        #
        # So dataProductType = cube
        #
        # PD 02-09-20
        # if those two files are images in the normal sense then it could make
        # sense to create separate planes with dataProductType = image that
        # end up with the correct (distinct) energy metadata.
        result = DataProductType.IMAGE
    else:
        mode = em.om.get('mode')
        obs_type = _get_obs_type(header)
        if mode is None:
            raise mc.CadcException(
                f'No mode information found for {em.om.get("filename")}'
            )
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
                    f'Mystery GPI mode {mode} for {em.om.get("filename")}'
                )
        elif (mode == 'imaging') or (
            obs_type is not None and obs_type == 'MASK'
        ):
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
    result = em.om.get('release')
    if result is not None and result.startswith('0001'):
        # because obs id GN-2008A-Q-39-69-015
        result = result.replace('0001', '2001')
    return result


def get_dec(header):
    """
    Get the declination. Rely on the JSON metadata, because it's all in
    the same units (degrees).

    :param header:  The FITS header for the current extension.
    :return: declination, or None if not found.
    """
    instrument = get_instrument()
    if instrument is em.Inst.HOKUPAA:
        ra, dec = _get_sky_coord(header, 'RA', 'DEC')
        result = dec
    elif instrument is em.Inst.OSCIR:
        ra, dec = _get_sky_coord(header, 'RA_TEL', 'DEC_TEL')
        result = dec
    elif instrument in [em.Inst.BHROS, em.Inst.TEXES]:
        # bHROS, TEXES ra/dec not in json
        # NIFS ra/dec not reliable in json
        # DB - 18-05-21
        # NIFS removed from 'if' statement
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
    instrument = get_instrument()
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


def get_meta_release(parameters):
    """
    Determine the metadata release date (Observation and Plane-level).

    :param parameters:  A dictionary container the FITS header for the
        current extension, as well as the URI for the .
    :return: The Observation/Plane release date, or None if not found.
    """
    uri = parameters.get('uri')
    if uri is None:
        raise mc.CadcException('uri missing from parameters.')

    # make sure the metadata is for the correct plane/file
    # combination - this location happens to be the first function called
    # during blueprint evaluation, which is why reset is
    # called here
    file_id = GemName.remove_extensions(mc.CaomName(uri).file_name)
    em.om.reset_index(file_id)

    header = parameters.get('header')
    if header is None:
        # GenericParser, so no headers retrieved from archive.gemini.edu,
        # probably a 403 being returned by the site, assume proprietary
        meta_release = em.om.get('release')
    else:
        # DB 21-08-19
        # If PROP_MD is T, use JSON ‘release’ value for metadata release date.
        # If no PROP_MD present or value is F use the JSON ut_datetime value.
        prop_md = header.get('PROP_MD')
        if prop_md is None or prop_md is False or prop_md == 'F':
            meta_release = em.om.get('ut_datetime')
        else:
            meta_release = em.om.get('release')
    return meta_release


def get_obs_intent(header):
    """
    Determine the Observation intent.

    :param header:  The FITS header for the current extension.
    :return: The Observation intent, or None if not found.
    """
    result = ObservationIntentType.CALIBRATION
    # DB 01-04-21
    # PINHOLE is CALIBRATION
    # DB 03-06-21
    # OBSCLASS = dayCal datasets should all have an intent of calibration.
    cal_values = [
        'GCALflat',
        'Bias',
        'BIAS',
        'Twilight',
        'Ar',
        'FLAT',
        'flat',
        'ARC',
        'Domeflat',
        'DARK',
        'dark',
        'gcal',
        'ZERO',
        'SLIT',
        'slit',
        'PINHOLE',
        'dayCal',
    ]
    dl = em.om.get('data_label')
    if dl is None and header is not None:
        dl = header.get('DATALAB')
    lookup = _get_obs_class(header)
    logging.debug(f'observation_class is {lookup} for {dl}')
    if lookup is None:
        type_lookup = _get_obs_type(header)
        logging.debug(f'observation_type is {type_lookup} for {dl}')
        if type_lookup is None:
            data_label = _get_data_label()
            logging.debug(f'data_label is {data_label}')
            if (
                data_label is None
                or (data_label is not None and '-CAL' not in data_label)
                and header is not None
            ):
                object_value = header.get('OBJECT')
                logging.debug(f'object_value is {object_value} for {dl}')
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
                instrument = get_instrument()
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
    instrument = get_instrument()
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
    elif instrument is em.Inst.FLAMINGOS:
        ignore, obs_type = _get_flamingos_mode(header)
        result = obs_type
    elif instrument is em.Inst.F2:
        # DB 03-06-21
        # check the 'types' JSON value
        types = em.om.get('types')
        if 'DARK' in types:
            result = 'DARK'
    return result


def get_proposal_id(header):
    """
    Determine the Proposal ID.

    :param header:  The FITS header for the current extension.
    :return: The proposal id from Gemini JSON metadata, or None if not found.
    """
    return em.om.get('program_id')


def get_provenance_keywords(uri):
    """
    DB https://github.com/opencadc-metadata-curation/gem2caom2/issues/12
    Currently there is no CAOM2 metadata that enables a user to distinguish
    different spectroscopic modes of GMOS-N/S data. e.g. long-slit vs. IFU
    vs. MOS (multi-object spectroscopy).

    To enable this with at least a TAP query it would be useful to modify the
    gem2caom2 code to add an Instrument.keywords value for GMOS-N/S spectra.

    The jsonsummary 'mode' value is likely sufficient for this. It is
    supposed to provide values of "imaging, spectroscopy, LS (Longslit
    Spectroscopy), MOS (Multi Object Spectroscopy) or IFS (Integral Field
    Spectroscopy)". There is likely no reason to change these values but
    simply use them for the value of Instrument.keywords.
    :param uri:
    :return:
    """
    return em.om.get('mode')


def get_provenance_last_executed(parameters):
    def breakout(comments):
        result = None
        temp = comments.split('\n')
        if len(temp) > 6 and 'HST' in temp[6]:
            # go from HST to UTC
            result = mc.make_time(temp[6])
        return result

    return _get_provenance_breakout(parameters, breakout)


def get_provenance_producer(parameters):
    result = None
    uri = parameters.get('uri')
    cal_level = get_calibration_level(uri)
    header = parameters.get('header')
    if cal_level in [
        CalibrationLevel.CALIBRATED,
        CalibrationLevel.PRODUCT,
        CalibrationLevel.ANALYSIS_PRODUCT,
    ]:
        comments = str(header.get('COMMENT'))
        result = comments.split('Processed by the')[1].split('|')[0]
    else:
        result = header.get('IMAGESWV')
    return result


def get_provenance_reference(parameters):
    def breakout(comments):
        return 'https://www.gemini.edu/instrumentation/graces/data-reduction'

    return _get_provenance_breakout(parameters, breakout)


def get_provenance_version(parameters):
    def breakout(comments):
        temp = comments.split('opera-')[1].split('build date')[0]
        return f'opera-{temp}'

    return _get_provenance_breakout(parameters, breakout)


def _get_provenance_breakout(parameters, fn):
    result = None
    uri = parameters.get('uri')
    cal_level = get_calibration_level(uri)
    if cal_level in [
        CalibrationLevel.CALIBRATED,
        CalibrationLevel.PRODUCT,
        CalibrationLevel.ANALYSIS_PRODUCT,
    ]:
        header = parameters.get('header')
        comments = str(header.get('COMMENT'))
        result = fn(comments)
    return result


def get_ra(header):
    """
    Get the right ascension. Rely on the JSON metadata, because it's all in
    the same units (degrees).

    :param header:  The FITS header for the current extension.
    :return: ra, or None if not found.
    """
    instrument = get_instrument()
    if instrument is em.Inst.HOKUPAA:
        ra, dec = _get_sky_coord(header, 'RA', 'DEC')
        result = ra
    elif instrument is em.Inst.OSCIR:
        ra, dec = _get_sky_coord(header, 'RA_TEL', 'DEC_TEL')
        result = ra
    elif instrument in [em.Inst.BHROS, em.Inst.TEXES]:
        # bHROS, TEXES: ra/dec not in json
        # DB - 05-03-19 - NIFS: ra/dec not reliable in json
        # DB - 18-05-21
        # NIFS removed from 'if' statement
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

    DB, 01-08-19
    Many calibration observations are acquired with the telescope parked and
    hence not tracking at sidereal rate.

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
    result = TargetType.FIELD
    spectroscopy = em.om.get('spectroscopy')
    instrument = get_instrument()
    if spectroscopy or instrument in [
        em.Inst.ALOPEKE,
        em.Inst.TEXES,
        em.Inst.ZORRO,
    ]:
        result = TargetType.OBJECT
    return result


def get_time_function_val(header):
    """
    Calculate the Chunk Time WCS function value, in 'mjd'.

    :param header:  The FITS header for the current extension (not used).
    :return: The Time WCS value from JSON Summary Metadata.
    """
    instrument = get_instrument()
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
    obs_class = em.om.get('observation_class')
    if obs_class is None and header is not None:
        obs_class = header.get('OBSCLASS')
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
        result = mc.to_float(pix_scale) / 3600.0
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
        raise mc.CadcException(
            f'No mode information found for {em.om.get("filename")}'
        )
    if obs_type is None:
        # DB - Also, since I’ve found FLAMINGOS spectra if OBJECT keyword or
        # json value contains ‘arc’ (any case) then it’s an ARC observation
        # type
        #
        # For FLAMINGOS since OBS_TYPE seems to always be set to ‘Object’
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


def _get_sky_coord(header, ra_key, dec_key):
    ra_hours = header.get(ra_key)
    dec_hours = header.get(dec_key)
    if (
        ra_hours is None
        or dec_hours is None
        or ra_hours == 'INDEF'
        or dec_hours == 'INDEF'
        or ra_hours == 'Unknown'
        or dec_hours == 'Unknown'
    ):
        ra_deg = None
        dec_deg = None
    else:
        result = SkyCoord(
            f'{ra_hours} {dec_hours}',
            unit=(units.hourangle, units.deg),
        )
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
    elif 'dark' in object_value or (
        view_pos is not None and 'dark' in view_pos
    ):
        result = 'DARK'
    elif 'gcal' in object_value:
        result = 'FLAT'
    elif 'arc' in object_value:
        result = 'ARC'
    elif 'slit' in object_value:
        # DB 22-02-21
        # These are images that show the slit location so I think it’s best
        # to add a new OBSTYPE of SLIT (sort of like MASK for GMOS-N/S).
        result = 'SLIT'
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
    if (
        'flat ' in object_value
        or 'dark ' in object_value
        or 'arc' in object_value
        or 'comp' in object_value
        or 'lamp' in object_value
        or 'comparison' in object_value
        or 'slit' in object_value
    ):
        return True
    return False


def _is_gmos_mask(header):
    result = False
    instrument = get_instrument()
    if instrument in [em.Inst.GMOSS, em.Inst.GMOSN, em.Inst.GMOS]:
        obs_type = _get_obs_type(header)
        if obs_type == 'MASK':
            result = True
    return result


def accumulate_fits_bp(bp, file_id, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model
    Observation level."""
    logging.debug(f'Begin accumulate_fits_bp for {file_id}.')
    em.get_obs_metadata(file_id)

    meta_producer = mc.get_version(APPLICATION)
    bp.set('Observation.type', 'get_obs_type(header)')
    bp.set('Observation.intent', 'get_obs_intent(header)')
    bp.set('Observation.metaProducer', meta_producer)
    bp.set('Observation.metaRelease', 'get_meta_release(parameters)')
    bp.set('Observation.target.type', 'get_target_type(uri)')
    bp.set('Observation.target.moving', 'get_target_moving(header)')
    bp.set('Observation.proposal.id', 'get_proposal_id(header)')

    bp.clear('Observation.algorithm.name')
    instrument = get_instrument()
    if instrument in [em.Inst.GMOSN, em.Inst.GMOSS, em.Inst.GMOS]:
        bp.set(
            'Observation.instrument.keywords', 'get_provenance_keywords(uri)'
        )
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
    bp.set('Plane.metaProducer', meta_producer)
    bp.set('Plane.metaRelease', 'get_meta_release(parameters)')
    bp.set('Plane.dataRelease', 'get_data_release(header)')

    bp.set('Plane.provenance.name', 'Gemini Observatory Data')
    bp.set('Plane.provenance.project', 'Gemini Archive')
    # Add IMAGESWV for GRACES
    bp.add_fits_attribute('Plane.provenance.producer', 'IMAGESWV')
    bp.set_default('Plane.provenance.producer', 'Gemini Observatory')
    if instrument in [em.Inst.ALOPEKE, em.Inst.ZORRO]:
        bp.set(
            'Plane.provenance.reference',
            f'http://archive.gemini.edu/searchform/filepre={file_id}.fits',
        )
    elif instrument is not em.Inst.TEXES:
        data_label = _get_data_label()
        bp.set(
            'Plane.provenance.reference',
            f'http://archive.gemini.edu/searchform/{data_label}',
        )

    if instrument is em.Inst.GRACES:
        bp.set(
            'Plane.provenance.lastExecuted',
            'get_provenance_last_executed(parameters)',
        )
        bp.set(
            'Plane.provenance.producer',
            'get_provenance_producer(parameters)',
        )
        bp.set(
            'Plane.provenance.reference',
            'get_provenance_reference(parameters)',
        )
        bp.set(
            'Plane.provenance.version',
            'get_provenance_version(parameters)',
        )

    bp.set('Artifact.metaProducer', meta_producer)
    bp.set('Artifact.productType', 'get_art_product_type(header)')
    bp.set('Artifact.contentChecksum', f'md5:{em.om.get("data_md5")}')
    bp.set('Artifact.contentLength', em.om.get('data_size'))
    bp.set('Artifact.contentType', 'application/fits')
    # always see the metadata, see the data only when it's public
    bp.set('Artifact.releaseType', 'data')
    bp.set('Artifact.uri', uri)

    if instrument is em.Inst.CIRPASS:
        bp.set_default('Observation.telescope.name', 'Gemini-South')
    mode = em.om.get('mode')
    if not (
        instrument
        in [
            em.Inst.GPI,
            em.Inst.PHOENIX,
            em.Inst.HOKUPAA,
            em.Inst.OSCIR,
            em.Inst.BHROS,
            em.Inst.TRECS,
        ]
        or (
            instrument is em.Inst.GRACES
            and mode is not None
            and mode != 'imaging'
        )
    ):
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

    bp.set('Chunk.metaProducer', meta_producer)
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
    if instrument in [em.Inst.ALOPEKE, em.Inst.ZORRO]:
        bp.clear('Chunk.time.axis.function.naxis')
        bp.add_fits_attribute('Chunk.time.axis.function.naxis', 'NAXIS3')
        bp.set_default('Chunk.time.axis.function.naxis', 1)

    bp.set('Chunk.time.axis.function.delta', 'get_time_delta(header)')
    bp.set('Chunk.time.axis.function.refCoord.pix', '0.5')
    bp.set(
        'Chunk.time.axis.function.refCoord.val',
        'get_time_function_val(header)',
    )

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
    else:
        current_product_id = None

    if headers is None:
        logging.info(
            f'Returning an un-modified observation '
            f'{observation.observation_id}.'
        )
        return observation

    if observation.instrument.name == 'oscir':
        # for these observations:
        # GN-2001A-C-16-3-016
        # GN-2001A-C-2-14-015
        # GN-2001A-C-2-2-002
        # GN-2001A-C-2-3-003
        # GN-2001A-C-2-4-004
        # GN-2001A-C-2-5-005
        # GN-2001A-C-2-6-006
        # GN-2001A-C-2-7-007
        # GN-2001A-C-2-8-009
        # GN-2001A-C-2-9-010
        observation.instrument = Instrument(name='OSCIR')
    instrument = em.Inst(observation.instrument.name)

    # processed files
    if cc.is_composite(headers) and not isinstance(
        observation, DerivedObservation
    ):
        observation = _update_composite(
            observation, instrument, current_product_id
        )

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
                f'{instrument}: no image data for '
                f'{observation.observation_id}. Cannot build an observation.'
            )
            return None

    config = mc.Config()
    config.get_executors()
    try:
        for plane in observation.planes.values():
            delete_list = []
            if (
                current_product_id is not None
                and current_product_id != plane.product_id
            ):
                continue

            for artifact in plane.artifacts.values():
                _should_artifact_be_deleted(artifact, config, delete_list)
                if GemName.is_preview(artifact.uri):
                    continue

                caom_name = mc.CaomName(artifact.uri)
                file_id = GemName.remove_extensions(
                    mc.CaomName(caom_name.uri).file_name
                )
                em.om.reset_index(file_id)
                processed = ofr.is_processed(caom_name.file_name)
                if instrument in [
                    em.Inst.MICHELLE,
                    em.Inst.TRECS,
                    em.Inst.GNIRS,
                ]:
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
                        artifact,
                        headers,
                        instrument,
                        observation.observation_id,
                    )

                for part in artifact.parts:

                    if part == '2' and instrument is em.Inst.GPI:
                        # GPI data sets have two extensions. First is science
                        # image (with WCS), second is data quality for each
                        # pixel (no WCS).
                        logging.info(
                            f'GPI: Setting chunks to None for part {part} '
                            f'for {observation.observation_id}'
                        )
                        artifact.parts[part].chunks = TypedList(
                            Chunk,
                        )
                        continue
                    for c in artifact.parts[part].chunks:
                        header = headers[int(part)]

                        # energy WCS
                        filter_name = get_filter_name(
                            headers[0],
                            header,
                            observation.observation_id,
                            instrument,
                        )
                        x = instruments.instrument_factory(
                            instrument, headers, header
                        )
                        x.filter_name = filter_name
                        if x.reset_energy(observation.type):
                            cc.reset_energy(c)
                        else:
                            x.chunk = c
                            x.data_product_type = plane.data_product_type
                            x.obs_id = observation.observation_id
                            x.update_energy()

                        # position WCS
                        mode = em.om.get('mode')
                        if x.reset_position(headers, observation.type):
                            logging.debug(
                                f'Setting Spatial WCS to None for '
                                f'{observation.observation_id}'
                            )
                            cc.reset_position(c)
                        else:
                            if instrument in [
                                em.Inst.PHOENIX,
                                em.Inst.HOKUPAA,
                                em.Inst.OSCIR,
                            ] or (
                                instrument is em.Inst.GRACES
                                and mode is not None
                                and mode != 'imaging'
                            ):
                                _update_chunk_position(
                                    c,
                                    header,
                                    instrument,
                                    int(part),
                                    observation.observation_id,
                                )
                            elif instrument in [
                                em.Inst.BHROS,
                                em.Inst.CIRPASS,
                                em.Inst.TEXES,
                            ]:
                                _update_chunk_position(
                                    c,
                                    headers[0],
                                    instrument,
                                    int(part),
                                    observation.observation_id,
                                )
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
                                    c,
                                    headers[0],
                                    instrument,
                                    int(part),
                                    observation.observation_id,
                                    n_axis1,
                                    n_axis2,
                                )
                            elif instrument is em.Inst.GPI:
                                _update_chunk_position(
                                    c,
                                    headers[1],
                                    instrument,
                                    int(part),
                                    observation.observation_id,
                                )
                                if part == '1':
                                    # equinox information only available from
                                    # 0th header
                                    equinox = headers[0].get('TRKEQUIN')
                                    if (
                                        equinox is not None
                                        and 1800.0 <= equinox <= 2500.0
                                    ):
                                        c.position.equinox = equinox
                                    else:
                                        # DB 07-06-21
                                        # No spatial WCS in these cases.
                                        cc.reset_position(c)

                            elif instrument is em.Inst.FLAMINGOS:
                                _update_chunk_position_flamingos(
                                    c, header, observation.observation_id
                                )
                            elif instrument is em.Inst.TRECS:
                                # DB 22-08-19
                                # For the file with CTYPE1 = 0 in HDU 1, all
                                # of the other WCS info appears to be fine
                                # (i.e. identical with the primary header).
                                # Use the HDU 0 values.
                                _update_chunk_position_trecs(
                                    c,
                                    headers,
                                    int(part),
                                    observation.observation_id,
                                )
                            elif instrument in [
                                em.Inst.ALOPEKE,
                                em.Inst.ZORRO,
                            ]:
                                _update_chunk_position_fox(
                                    c, observation.observation_id
                                )
                            elif instrument is em.Inst.NIRI:
                                _update_chunk_position_niri(
                                    c,
                                    headers,
                                    observation.observation_id,
                                    int(part),
                                )

                        # time WCS
                        if instrument is em.Inst.F2:
                            _update_chunk_time_f2(
                                c, observation.observation_id
                            )

                        # DB - 05-06-20
                        # That’s a composite observation (but with no way of
                        # determining the 4 members from the header) and it is
                        # an extracted spectrum and (despite some header info
                        # suggesting otherwise) has no wavelength scale.
                        # LINEAR must be a reference to the fact that the
                        # spacing of each pixel is constant.  Since the scale
                        # is simply in pixels….  CTYPE1 will refer to the
                        # energy axis.   Would likely make more sense for a
                        # value of PIXEL.
                        if instrument is em.Inst.PHOENIX:
                            ctype = headers[0].get('CTYPE1')
                            if ctype is not None and ctype in [
                                'LINEAR',
                                'PIXEL',
                            ]:
                                c.naxis = None
                                c.position_axis_1 = None
                                c.position_axis_2 = None
                                c.time_axis = None
                                c.energy_axis = None

                        if instrument in [em.Inst.ALOPEKE, em.Inst.ZORRO]:
                            _update_chunk_time_fox(
                                c, header, observation.observation_id
                            )
                        if c.naxis is not None and c.naxis <= 2:
                            if c.position_axis_1 is None:
                                c.naxis = None
                            c.time_axis = None

                        # DB - 01-02-21
                        # OSCIR data should all have NAXIS=6 in the header.
                        # Axis 1/2 are position axes.  Axis 3 is the number of
                        # chop positions (1 or usually = 2), axis 5 is the
                        # number of nod positions (1 or usually = 2), axis
                        # 4 gives the number of ‘savesets’ per nod position
                        # and axis 6 gives the number of ‘nod sets’.
                        #
                        # For that particular example you gave me I think the
                        # NAXIS6 = 1 is incorrect.  Think it should equal the
                        # value of NODSETS.  And TOTFRMS = 360 should equal
                        # naxis1 x naxis2 x naxis3 x naxis4.
                        #
                        # Spatial cutouts would use axes 1 and 2, although
                        # the precise positions of each chop/nod position are
                        # not captured by the CAOM2 data.
                        if (
                            instrument is em.Inst.OSCIR
                            and c.naxis is not None
                            and c.naxis == 6
                        ):
                            c.naxis = 2
                            c.time_axis = None
                        if (
                            instrument in [em.Inst.MICHELLE, em.Inst.TRECS]
                            and c.naxis == 4
                        ):
                            if c.position is None:
                                c.naxis = None
                                if c.time_axis is not None:
                                    c.time_axis = None
                            else:
                                if c.naxis == 4:
                                    c.naxis = 2
                                    c.time_axis = None

                        # DB 04-17-21
                        # BIASes and DARKs for all instruments should ignore
                        # spatial and spectral wcs - clean up associated
                        # axes
                        # also fix a very specific edge case where cal files
                        # have useless WCS information for the purposes of
                        # CAOM2.4 axis checks, and the corresponding
                        # cutouts. No spatial wcs means invalid chunk.naxis
                        # value, so set that to None, which then
                        # invalidates the chunk.time_axis value
                        # DB 06-01-21 - no calibration file cutouts to
                        # support, so removing this axis information is
                        # not removing downstream functionality
                        if (
                            c.naxis == 3
                            and c.position is None
                            and c.time is not None
                        ):
                            if c.time.axis.function.naxis == 1:
                                c.naxis = None
                                c.time_axis = None
                            else:
                                c.naxis = 1
                                c.time_axis = 1

                if isinstance(observation, DerivedObservation):
                    values = cc.find_keywords_in_headers(
                        headers[1:], ['IMCMB']
                    )
                    repaired_values = _remove_processing_detritus(
                        values, observation.observation_id
                    )
                    cc.update_plane_provenance_from_values(
                        plane,
                        _repair_provenance_value,
                        repaired_values,
                        COLLECTION,
                        observation.observation_id,
                    )

                if (
                    processed
                    or isinstance(observation, DerivedObservation)
                    or instrument is em.Inst.TEXES
                ) and 'jpg' not in caom_name.file_name:
                    # not the preview artifact
                    if plane.provenance is not None:
                        if instrument is not em.Inst.GRACES:
                            plane.provenance.reference = (
                                f'http://archive.gemini.edu/searchform/'
                                f'filepre={caom_name.file_name}'
                            )

            program = em.get_pi_metadata(observation.proposal.id)
            if program is not None:
                observation.proposal.pi_name = program['pi_name']
                observation.proposal.title = program['title']

            temp = list(set(delete_list))
            for entry in temp:
                logging.warning(
                    f'Removing artifact {entry} from observation '
                    f'{observation.observation_id}, plane {plane.product_id}.'
                )
                plane.artifacts.pop(entry)

        if isinstance(observation, DerivedObservation):
            cc.update_observation_members(observation)

        em.value_repair.repair(observation)
    except Exception as e:
        logging.error(
            f'Error {e} for {observation.observation_id} instrument '
            f'{instrument}'
        )
        tb = traceback.format_exc()
        logging.debug(tb)
        raise mc.CadcException(e)
    logging.debug('Done update.')
    return observation


def _should_artifact_be_deleted(artifact, config, delete_list):
    if config.features.supports_latest_client:
        if artifact.uri.startswith('gemini'):
            if 'GEMINI' not in artifact.uri:
                delete_list.append(artifact.uri)
        if artifact.uri.startswith('ad'):
            delete_list.append(artifact.uri)


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
    if filter_name is None or (
        filter_name is not None and len(filter_name.strip()) == 0
    ):
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
        f'Filter names are {filter_name} for instrument {instrument} in '
        f'{obs_id}'
    )
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
    if ('AZEL_TARGET' in types and ra is None) or (
        instrument is em.Inst.GNIRS and ra_json is None and dec_json is None
    ):
        # DB - 02-04-19 - Az-El coordinate frame likely means the telescope
        # was parked or at least not tracking so spatial information is
        # irrelevant.

        # DB - 09-04-19 - AZEL_TARGET should likely be checked for all
        # datasets, and means the spatial WCS should be ignored. since this
        # generally means the telescope is not tracking and so spatial WCS
        # info isn’t relevant since the position is changing with time.

        # DB 24-04-19
        # Ignore spatial WCS if RA/Dec are not available for GNIRS.

        logging.info(f'{instrument}: Spatial WCS is None for {obs_id}')
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


def _update_chunk_position_fox(chunk, obs_id):
    logging.debug(f'Begin _update_chunk_position_fox for {obs_id}')
    if (
        chunk is not None
        and chunk.position is not None
        and chunk.position.axis is not None
        and chunk.position.axis.axis1.ctype == 'RA--TAN'
    ):
        chunk.position.axis.axis1.ctype = 'RA---TAN'
    logging.debug(f'End _update_chunk_position_fox.')


def _update_chunk_position_flamingos(chunk, header, obs_id):
    # DB - I see nothing in astropy that will do a transformation from crota
    # form to CD matrix, but this is it:

    # cd1_1 = cdelt1 * cos (crota1)
    # cd1_2 = -cdelt2 * sin (crota1)
    # cd2_1 = cdelt1 * sin (crota1)
    # cd2_2 = cdelt2 * cos (crota1)

    # Note that there is not a crota2 keyword (it would have the same value
    # as crota1 if it existed)
    if (
        chunk is not None
        and chunk.position is not None
        and chunk.position.axis is not None
        and chunk.position.axis.function is not None
    ):
        crval1 = header.get('CRVAL1')
        crval2 = header.get('CRVAL2')
        if 0.0 <= crval1 <= 360.0 and -90.0 <= crval2 <= 90.0:
            c_delt1 = header.get('CDELT1')
            c_delt2 = header.get('CDELT2')
            c_rota1 = header.get('CROTA1')
            if (
                c_delt1 is not None
                and c_delt2 is not None
                and c_rota1 is not None
            ):
                chunk.position.axis.function.cd11 = c_delt1 * math.cos(c_rota1)
                chunk.position.axis.function.cd12 = -c_delt2 * math.sin(
                    c_rota1
                )
                chunk.position.axis.function.cd21 = c_delt1 * math.sin(c_rota1)
                chunk.position.axis.function.cd22 = c_delt2 * math.cos(c_rota1)
            else:
                logging.info(
                    f'FLAMINGOS: Missing spatial wcs inputs for {obs_id}'
                )
                chunk.position.axis.function.cd11 = None
                chunk.position.axis.function.cd12 = None
                chunk.position.axis.function.cd21 = None
                chunk.position.axis.function.cd22 = None
        else:
            # DB 04-12-19
            # FLAMINGOS GS-CAL20020623-14-0080 02jun23.0080.fits
            # The header has “CRVAL1  =           3581.13808 ” which is
            # supposed to be in degrees and shouldn’t be > 360. Skip spatial
            # WCS if errors like this occur.
            logging.warning(
                f'FLAMINGOS: Spatial WCS set to None for {obs_id} '
                f'because CRVAL1 == {crval1} and CRVAL2 == {crval2}.'
            )
            cc.reset_position(chunk)
    else:
        logging.info(f'FLAMINGOS: Missing spatial wcs for {obs_id}')


def _update_chunk_position_niri(chunk, headers, obs_id, extension):
    logging.info(f'Begin _update_chunk_niri for {obs_id}')
    # DB 07-06-21
    # The extension CD values that are very, very close to 0 cause the
    # problems with Spatial WCS:
    # ERROR: spherepoly_from_array: a line segment overlaps or polygon too
    #        large
    # Try to use the primary values if this error occurs - there's an extra
    # '5' in the exponent
    if len(headers) > 1:
        pdu = headers[0]
        hdu0 = headers[1]
        pdu_cd1_1 = pdu.get('CD1_1')
        hdu0_cd1_1 = hdu0.get('CD1_1')
        if not math.isclose(pdu_cd1_1, hdu0_cd1_1) and math.isclose(
            pdu_cd1_1 * 1e-50, hdu0_cd1_1
        ):
            pdu['NAXIS1'] = hdu0.get('NAXIS1')
            pdu['NAXIS2'] = hdu0.get('NAXIS2')
            wcs_parser = WcsParser(pdu, obs_id, extension)
            if chunk is None:
                chunk = Chunk()
            wcs_parser.augment_position(chunk)
            if chunk.position is not None:
                chunk.position_axis_1 = 1
                chunk.position_axis_2 = 2
                chunk.position.coordsys = pdu.get('FRAME')
                chunk.position.equinox = mc.to_float(pdu.get('EQUINOX'))
    logging.info('End _update_chunk_niri')


def _update_chunk_position_trecs(chunk, headers, extension, obs_id):
    if len(headers) > 1:
        ctype1 = headers[extension].get('CTYPE1')
        if isinstance(ctype1, str):
            # value repair for a small subset of TReCS files  :(
            # test is rS20060306S0090, GS-2005B-Q-10-63-003
            if ctype1 == '0':
                headers[extension]['CTYPE1'] = 'RA---TAN'
            wcs_parser = WcsParser(headers[extension], obs_id, extension)
        else:
            wcs_parser = WcsParser(headers[0], obs_id, extension)
    else:
        wcs_parser = WcsParser(headers[0], obs_id, extension)

    if chunk is None:
        chunk = Chunk()

    wcs_parser.augment_position(chunk)
    if chunk.position is not None:
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
        chunk.position.coordsys = headers[0].get('FRAME')
        chunk.position.equinox = mc.to_float(headers[0].get('EQUINOX'))


def _update_chunk_position(
    chunk, header, instrument, extension, obs_id, n_axis1=None, n_axis2=None
):
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

    # DB - 21-08-19
    # GN-2001A-DD-2-1-537:  Much of the header info seems to be incorrect so
    # no WCS.  But what’s odd is that searching Gemini’s archive for this data
    # label returns 3 files with different file names:  23, 24 and 25
    # February.  EQUINOX only seems to be incorrect in 01FEB23_537.fits.
    #
    # GN-2001A-DD-2-3-1060:  header has UNKNOWN for values of coordinates so
    # no spatial WCS

    if instrument is em.Inst.HOKUPAA:
        equinox = mc.to_float(header.get('EQUINOX'))
        if not 1800.0 <= equinox <= 2500.0:
            logging.warning(
                f'EQUINOX value is wrong ({equinox}), no spatial WCS for '
                f'{obs_id}'
            )
            return

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
        # DB 18-05-21
        # EQUINOX and RADESYS values should be fine for all NIFS files since
        # that’s what is used for all other Gemini data as far as I’m aware.
        header['EQUINOX'] = 2000.0
        header['RADESYS'] = 'FK5'
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
        chunk.position.equinox = mc.to_float(header.get('EQUINOX'))
    elif instrument is em.Inst.BHROS:
        chunk.position.coordsys = header.get('TRKFRAME')
    elif instrument is em.Inst.GPI:
        chunk.position.coordsys = header.get('RADESYS')

    logging.debug('End _update_chunk_position')


def _update_chunk_time_f2(chunk, obs_id):
    """F2 FITS files have a CD3_3 element that's not supported by fits2caom2,
    so using the blueprint will not work to adjust that value. Set delta
    specifically here."""
    logging.debug(f'Begin _update_chunk_time_f2 {obs_id}')
    mc.check_param(chunk, Chunk)
    if (
        chunk.time is not None
        and chunk.time.axis is not None
        and chunk.time.axis.function is not None
    ):
        chunk.time.axis.function.delta = get_time_delta(None)
        logging.info(f'F2: Updated time delta for {obs_id}')
    logging.debug(f'End _update_chunk_time_f2 {obs_id}')


def _update_chunk_time_fox(chunk, header, obs_id):
    """
    DB 02-09-20
    Exposure time using JSON values isn’t correct.  I know that for this
    example Gemini shows the exposure time is 0.02 seconds but there are
    1000 x 0.02-second exposures in the cube.  The keyword EXPOSURE gives the
    total exposure time (in seconds), time.exposure, or 20 in this case while
    the json exposure_time should be the time.resolution.
    """
    logging.debug(f'Begin _update_chunk_time_fox {obs_id}')
    mc.check_param(chunk, Chunk)
    if chunk.time is not None:
        chunk.time.exposure = header.get('EXPOSURE')
        # chunk.time.resolution already set by blueprint
    logging.debug(f'End _update_chunk_time_fox {obs_id}')


def _update_chunk_time_gmos(chunk, obs_id):
    """"""
    logging.debug(f'Begin _update_chunk_time_gmos {obs_id}')
    mc.check_param(chunk, Chunk)
    if chunk.time is not None and chunk.time.axis is not None:
        mjd_start = get_time_function_val(None)
        start = RefCoord(0.5, mjd_start.value)
        end = RefCoord(1.5, mjd_start.value)
        chunk.time.exposure = 0.0
        chunk.time.axis.range = CoordRange1D(start, end)
        chunk.time.axis.function = None
    logging.debug(f'End _update_chunk_time_gmos {obs_id}')


def _update_composite(obs, instrument, current_product_id):
    if instrument is em.Inst.TRECS:
        if current_product_id is not None and (
            current_product_id.startswith('rS')
            or current_product_id.startswith('rN')
        ):
            # DB 02-06-20
            # processed TReCS files in Gemini's archive are derived by
            # combining the NNODSETS x NSAVSETS contained within a single
            # unprocessed image into a simpler image array.
            #
            # SGo - this means ignoring the IMCMB keywords that are an
            # artifact of that, which is how Composite construction is
            # otherwise determined.
            result = obs
    else:
        result = cc.change_to_composite(obs)
        logging.info(f'{obs.observation_id} is a Composite Observation.')
    return result


def _remove_processing_detritus(values, obs_id):
    """
    There are several naming patterns in the provenance for processed files.
    Try to extract meaningful raw file names.

    Remove duplicates, so values used for provenance in multiple forms are
    only looked up once.
    :param values: A list of IMCMB* keyword values.
    :param obs_id: str for logging information
    :return: A list of unique provenance file names, such as may be found
        at Gemini. The list may be empty.
    """
    logging.debug(f'Begin _remove_processing_detritus for {obs_id}')
    result = []
    for value in values:
        # e.g.
        # IMCMB001 = 'tmpimgwsk9476kd_5.fits[SCI,1]'
        # tmpfile22889S20141226S0203.fits[SCI,1]
        # IMCMB001= 'rawdir$2004may20_0048.fits'

        if 'N' in value:
            temp = 'N' + value.split('N', 1)[1]
        elif 'S' in value:
            x = value.split('S')
            if len(x) == 2 and '[SCI' in value:
                logging.warning(f'Unrecognized IMCMB value {value}')
                continue
            else:
                temp = 'S' + value.split('S', 1)[1]
        elif '$' in value:
            temp = value.split('$', 1)[1]
        else:
            logging.warning(f'Unrecognized IMCMB value {value}')
            continue

        if '_' in temp and (temp.startswith('S') or temp.startswith('N')):
            temp1 = temp.split('_')[0]
        elif '.fits' in temp:
            temp1 = temp.split('.fits')[0]
        elif '[SCI' in temp:
            temp1 = temp.split('[SCI')[0]
        else:
            logging.warning(f'Failure to repair {temp}')
            continue

        result.append(temp1[:14])

    logging.debug('End _remove_processing_detritus')
    return list(set(result))


def _repair_provenance_value(value, obs_id):
    logging.debug(f'Being _repair_provenance_value for {obs_id}.')
    prov_file_id = value
    try:
        prov_obs_id = em.get_gofr().get_obs_id(prov_file_id)
    except mc.CadcException as e:
        # the file id probably does not exist at Gemini, ignore, because
        # it's provenance
        logging.warning(f'Failed to find {prov_file_id} at archive.gemini.edu')
        # DB 01-06-21 - use not found for the DATALAB/observationID value
        # so it's easy to find in the database and let Gemini know.
        prov_obs_id = 'not_found'
    logging.debug(
        f'End _repair_provenance_value. {prov_obs_id} {prov_file_id}'
    )
    return prov_obs_id, prov_file_id


def _build_blueprints(uris):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param uris The list of artifact URIs for the files to be processed.
    """
    module = importlib.import_module(__name__)
    blueprints = {}
    for uri in uris:
        blueprint = ObsBlueprint(module=module)
        if not GemName.is_preview(uri):
            file_id = GemName.remove_extensions(mc.CaomName(uri).file_name)
            accumulate_fits_bp(blueprint, file_id, uri)
        blueprints[uri] = blueprint
    return blueprints


def _get_uris(args):
    result = []
    if args.lineage:
        for ii in args.lineage:
            ignore, temp = mc.decompose_lineage(ii)
            result.append(temp)
    elif args.local:
        config = mc.Config()
        config.get_executors()
        name_builder = GemObsIDBuilder(config)
        for ii in args.local:
            result.append(name_builder.build(os.path.basename(ii)).file_uri)
    else:
        raise mc.CadcException(f'Could not define uri from these args {args}')
    return result


def to_caom2():
    try:
        args = get_gen_proc_arg_parser().parse_args()
        uris = _get_uris(args)
        blueprints = _build_blueprints(uris)
        result = gen_proc(args, blueprints)
        return result
    except Exception as e:
        logging.error(traceback.format_exc())
        raise e


def gem_main_app():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        result = to_caom2()
        logging.debug(f'Done {APPLICATION} processing.')
        sys.exit(result)
    except Exception as e:
        logging.error(f'Failed {APPLICATION} execution for {args}.')
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)
