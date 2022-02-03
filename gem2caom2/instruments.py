# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
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
#  : 4 $
#
# ***********************************************************************
#

import logging
import math
import re

from caom2 import DataProductType, CoordAxis1D, RefCoord, Axis, CoordRange1D
from caom2 import SpectralWCS, Chunk, ObservationIntentType, TargetType
from caom2 import ProductType, CalibrationLevel
from caom2utils import fits2caom2
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from gem2caom2 import svofps
from gem2caom2.obs_metadata import json_lookup
from gem2caom2.util import Inst


__all__ = ['instrument_factory', 'InstrumentType']


# DB 01-04-21
# PINHOLE is CALIBRATION
# DB 03-06-21
# OBSCLASS = dayCal datasets should all have an intent of calibration.
CAL_VALUES = [
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


# DB - 04-02-19 - strip out anything with 'pupil' as it doesn't affect
# energy transmission
#
# DB 24-04-19
# Other instruments occasionally have ND filters in the beam.
FILTER_VALUES_TO_IGNORE = ['open', 'invalid', 'pupil', 'clear', 'nd']


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
# IGRINS - DB - 21-07-21 - hard-code 5" FOV
RADIUS_LOOKUP = {
    Inst.GPI: 2.8 / 3600.0,  # units are arcseconds
    Inst.GRACES: 1.2 / 3600.0,
    Inst.PHOENIX: 5.0 / 3600.0,
    Inst.OSCIR: 0.0890 / 3600.0,
    Inst.HOKUPAA: 4.0 / 3600.0,
    Inst.BHROS: 0.9 / 3600.0,
    Inst.NIFS: 3.0 / 3600.0,
    Inst.TEXES: 5.0 / 3600.0,
    Inst.IGRINS: 5.0 / 3600.0,
}
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


class InstrumentType:
    def __init__(self, name):
        self._name = name
        self._chunk = None
        self._data_label = None
        self._data_product_type = None
        self._fm = None
        self._filter_name = None
        self._extension = None
        self._headers = None
        self._mode = None
        self._obs_class = None
        self._obs_id = None
        self._obs_type = None
        self._svo_session = None
        self._logger = logging.getLogger(name.value)

    def _another_init_that_i_need_to_figure_out_where_to_call(self, header):
        self._get_data_label(header)
        self._get_obs_type(header)
        self._get_obs_class(header)

    @property
    def chunk(self):
        return self._chunk

    @chunk.setter
    def chunk(self, value):
        self._chunk = value

    @property
    def data_product_type(self):
        return self._data_product_type

    @data_product_type.setter
    def data_product_type(self, value):
        self._data_product_type = value

    @property
    def extension(self):
        return self._extension

    @extension.setter
    def extension(self, value):
        self._extension = value

    @property
    def filter_name(self):
        return self._filter_name

    @filter_name.setter
    def filter_name(self, value):
        self._filter_name = value
        if value == 'NONE':
            self._filter_name = None

    @property
    def fm(self):
        return self._fm

    @fm.setter
    def fm(self, value):
        self._fm = value

    @property
    def headers(self):
        return self._headers

    @headers.setter
    def headers(self, value):
        self._headers = value

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, value):
        self._mode = value

    @property
    def obs_id(self):
        return self._obs_id

    @obs_id.setter
    def obs_id(self, value):
        self._obs_id = value

    def _get_data_label(self, header):
        self._data_label = json_lookup.get('data_label')
        if self._data_label is None and header is not None:
            self._data_label = header.get('DATALAB')

    def _get_obs_class(self, header):
        """Common location to lookup OBSCLASS from the FITS headers, and if
        it's not present, to lookup observation_class from JSON summary
        metadata."""
        self._obs_class = json_lookup.get('observation_class')
        if self._obs_class is None and header is not None:
            self._obs_class = header.get('OBSCLASS')

    def _get_obs_type(self, header):
        """Common location to lookup OBSTYPE from the FITS headers, and if
        it's not present, to lookup observation_type from JSON summary
        metadata."""
        self._obs_type = json_lookup.get('observation_type')
        if self._obs_type is None:
            self._obs_type = header.get('OBSTYPE')

    def build_chunk_energy(self):
        # If n_axis=1 (as I guess it will be for all but processes GRACES
        # spectra now?) that means crpix=0.5 and the corresponding crval would
        # central_wl - bandpass/2.0 (i.e. the minimum wavelength).   It is fine
        # if you instead change crpix to 1.0.   I guess since the ‘range’ of
        # one pixel is 0.5 to 1.5.

        axis = CoordAxis1D(axis=Axis(ctype='WAVE', cunit='um'))
        ref_coord1 = RefCoord(0.5, self.fm.central_wl - self.fm.bandpass / 2.0)
        ref_coord2 = RefCoord(1.5, self.fm.central_wl + self.fm.bandpass / 2.0)
        axis.range = CoordRange1D(ref_coord1, ref_coord2)

        # DB - 14-02-19 value looks clearer (as two filters) with a space on
        # both sides of the ‘+’.
        bandpass_name = (
            None if (
                    self.filter_name is None or len(self.filter_name) == 0
            ) else self.filter_name.replace('+', ' + ')
        )

        if bandpass_name is not None and 'empty' in bandpass_name:
            # DB 02-12-19 - But some files have filters open1-6 and empty_01.
            # So likely less confusing to remove ‘empty*’ completely.
            #
            # Remove the 'empty' string here, now that min/max wavelength
            # calculations have been completed.
            bandpass_name = bandpass_name.replace('empty', '')
            if len(bandpass_name) == 0:
                bandpass_name = None
            elif '+' in bandpass_name:
                bandpass_name = bandpass_name.replace('+', '').strip()

        if math.isclose(self.fm.central_wl, 0.0):
            energy = SpectralWCS(
                axis=CoordAxis1D(axis=Axis(ctype='WAVE', cunit='um')),
                specsys='TOPOCENT',
                ssyssrc=None,
                ssysobs=None,
                bandpass_name=bandpass_name,
            )
        else:
            energy = SpectralWCS(
                axis=axis,
                specsys='TOPOCENT',
                ssyssrc='TOPOCENT',
                ssysobs='TOPOCENT',
                bandpass_name=bandpass_name,
                resolving_power=self.fm.resolving_power,
            )
        self.chunk.energy = energy
        # no chunk energy is derived from FITS file axis metadata, so no
        # cutouts to support
        self.chunk.energy_axis = None

    def get_art_product_type(self, header):
        """
        Calculate the Artifact ProductType.

        If obsclass is unknown then CAOM2 ProductType is set to CALIBRATION.
        This should only effect early data when OBSCLASS was not in the JSON
        summary metadata or FITS headers.

        :param header:  The FITS header for the current extension.
        :return: The Artifact ProductType, or None if not found.
        """
        self._get_data_label(header)
        self._get_obs_type(header)
        self._get_obs_class(header)

        self._logger.debug(
            f'obs type is {self._obs_type} obs class is {self._obs_class} '
            f'for {self._data_label}'
        )
        if self._obs_type is not None and self._obs_type == 'MASK':
            result = ProductType.AUXILIARY
        elif self._obs_class is None:
            # 'unknown' is the value of observation_class for CIRPASS
            if (
                self._obs_type is not None and
                self._obs_type in ['OBJECT', 'unknown']
            ):
                if self._data_label is not None and 'CAL' in self._data_label:
                    result = ProductType.CALIBRATION
                else:
                    result = ProductType.SCIENCE
            else:
                result = ProductType.CALIBRATION
        elif self._obs_class is not None and self._obs_class == 'science':
            result = ProductType.SCIENCE
        else:
            result = ProductType.CALIBRATION
        return result

    def get_calibration_level(self, uri):
        result = CalibrationLevel.RAW_STANDARD
        reduction = json_lookup.get('reduction')
        if (
            reduction is not None and (
                'PROCESSED' in reduction or 'PREPARED' in reduction
            )
        ):
            result = CalibrationLevel.CALIBRATED
        return result

    def get_cd11(self, header, keyword='CDELT1'):
        result = header.get(keyword)
        if result is None:
            result = RADIUS_LOOKUP.get(self._name)
        return result

    def get_cd22(self, header):
        return self.get_cd11(header, 'CDELT2')

    def get_data_product_type(self, header):
        """
        Calculate the Plane DataProductType.

        :param header:  The FITS header for the current extension.
        :return: The Plane DataProductType, or None if not found.
        """
        mode = json_lookup.get('mode')
        self._get_obs_type(header)
        if mode is None:
            raise mc.CadcException(
                f'No mode information found for '
                f'{json_lookup.get("filename")}'
            )
        elif (mode == 'imaging') or (
                self._obs_type is not None and self._obs_type == 'MASK'
        ):
            result = DataProductType.IMAGE
        else:
            result = DataProductType.SPECTRUM
        return result

    def get_dec(self, header):
        """
        Get the declination. Rely on the JSON metadata, because it's all in
        the same units (degrees).

        :param header:  The FITS header for the current extension.
        :return: declination, or None if not found.
        """
        return json_lookup.get('dec')

    def get_exposure(self, header):
        """
        Calculate the exposure time. EXPTIME in the header is not always the
        total exposure time for some of the IR instruments.  For these EXPTIME
        is usually the individual exposure time for each chop/nod and a bunch
        of chops/nods are executed for one observation.  Gemini correctly
         allows for this in the json ‘exposure_time’ value.

        :param header:  The FITS header for the current extension (unused).
        :return: The exposure time, or None if not found.
        """
        return json_lookup.get('exposure_time')

    def get_filter_name(self):
        """
        Create the filter names for use by update_energy methods.

        :return: The filter names, or None if none found.
        """
        filter_name = json_lookup.get('filter_name')

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
        if filter_name is None or len(filter_name.strip()) == 0:
            result = InstrumentType._search_through_keys(
                self._headers[0], ['FILTER'], FILTER_VALUES_TO_IGNORE)
            filter_name = result
        self._logger.info(
            f'Filter names are "{filter_name}" in {self._data_label}'
        )
        self.filter_name = filter_name

    @staticmethod
    def _search_through_keys(header, search_keys, values_to_ignore):
        result = []
        for key in header.keys():
            for lookup in search_keys:
                if lookup in key:
                    value = header.get(key).lower()
                    ignore = False
                    for ii in values_to_ignore:
                        if ii.startswith(value) or value.startswith(ii):
                            ignore = True
                            break
                    if ignore:
                        continue
                    else:
                        result.append(header.get(key).strip())
        return '+'.join(result)

    def get_obs_intent(self, header):
        """
        Determine the Observation intent.

        :param header:  The FITS header for the current extension.
        :return: The Observation intent, or None if not found.
        """
        result = ObservationIntentType.CALIBRATION
        self._get_data_label(header)
        self._get_obs_class(header)
        self._logger.debug(
            f'observation_class is {self._obs_class} for {self._data_label}'
        )
        if self._obs_class is None:
            self._get_obs_type(header)
            self._logger.debug(
                f'observation_type is {self._obs_type} for {self._data_label}'
            )
            if self._obs_type is None:
                self._logger.debug(f'data_label is {self._data_label}')
                if (
                        self._data_label is None or (
                            self._data_label is not None and
                            '-CAL' not in self._data_label
                        )
                        and header is not None
                ):
                    object_value = header.get('OBJECT')
                    self._logger.debug(
                        f'object_value is {object_value} for '
                        f'{self._data_label}'
                    )
                    if object_value is not None:
                        if object_value in CAL_VALUES:
                            result = ObservationIntentType.CALIBRATION
                        else:
                            # check that an individual cal value is not part
                            # of the object_value
                            cal_value_found = False
                            for ii in CAL_VALUES:
                                if ii in object_value:
                                    result = ObservationIntentType.CALIBRATION
                                    cal_value_found = True
                                    break
                            if not cal_value_found:
                                result = ObservationIntentType.SCIENCE
                else:
                    if '-CAL' in self._data_label:
                        result = ObservationIntentType.CALIBRATION
                    else:
                        result = ObservationIntentType.SCIENCE
            else:
                if self._obs_type in CAL_VALUES:
                    result = ObservationIntentType.CALIBRATION
                else:
                    result = ObservationIntentType.SCIENCE
        elif 'science' in self._obs_class:
            result = ObservationIntentType.SCIENCE
        return result

    def get_obs_type(self, header):
        """
        Determine the Observation type.

        :param header:  The FITS header for the current extension.
        :return: The Observation type, or None if not found.
        """
        self._get_obs_type(header)
        result = self._obs_type
        self._get_obs_class(header)
        if (
            self._obs_class is not None and
            (self._obs_class == 'acq' or self._obs_class == 'acqCal')
        ):
            result = 'ACQUISITION'
        return result

    def get_ra(self, header):
        """
        Get the right ascension. Rely on the JSON metadata, because it's all in
        the same units (degrees).

        :param header:  The FITS header for the current extension.
        :return: ra, or None if not found.
        """
        return json_lookup.get('ra')

    def get_target_type(self):
        """
        Calculate the Target TargetType
        """
        result = TargetType.FIELD
        spectroscopy = json_lookup.get('spectroscopy')
        if spectroscopy:
            result = TargetType.OBJECT
        return result

    def get_time_delta(self, header):
        """
        Calculate the Time WCS delta.

        :param header: The FITS header for the current extension.
        :return: The Time delta, or None if none found.
        """
        exptime = self.get_exposure(header)
        if exptime is None:
            return None
        return mc.to_float(exptime) / (24.0 * 3600.0)

    def get_time_function_val(self, header):
        """
        Calculate the Chunk Time WCS function value, in 'mjd'.

        :param header:  The FITS header for the current extension (not used).
        :return: The Time WCS value from JSON Summary Metadata.
        """
        time_string = json_lookup.get('ut_datetime')
        return ac.get_datetime(time_string)

    def make_axes_consistent(self):
        # DB 04-17-21
        # BIASes and DARKs for all instruments should ignore spatial and
        # spectral wcs - clean up associated axes
        #
        # also fix a very specific edge case where cal files have useless WCS
        # information for the purposes of CAOM2.4 axis checks, and the
        # corresponding cutouts. No spatial wcs means invalid chunk.naxis
        # value, so set that to None, which then invalidates the
        # chunk.time_axis value DB 06-01-21 - no calibration file cutouts to
        # support, so removing this axis information is not removing
        # downstream functionality
        if (
            self.chunk.naxis == 3
            and self.chunk.position is None
            and self.chunk.time is not None
        ):
            if self.chunk.time.axis.function.naxis == 1:
                self.chunk.naxis = None
                self.chunk.time_axis = None
            else:
                self.chunk.naxis = 1
                self.chunk.time_axis = 1

        if (
            (self.chunk.naxis is not None and self.chunk.naxis <= 2) and not
            # the following exempts the Fox use case
            (self.chunk.naxis == 1 and self.chunk.time_axis == 1)
        ):
            if self.chunk.position_axis_1 is None:
                self.chunk.naxis = None
            self.chunk.time_axis = None

    def multiple_filter_lookup(self, lookup, wl_max=None):
        w_max = 10.0 if wl_max is None else wl_max
        w_min = 0.0
        for ii in self.filter_name.split('+'):
            if ii in lookup:
                wl_max = lookup[ii][2]
                wl_min = lookup[ii][1]
            else:
                msg = (
                    f'Unprepared for filter {ii} from {self.obs_id}'
                )
                if self._name == Inst.MICHELLE and (
                        ii.startswith('I') or (ii == 'Grid_T')
                ):
                    self._logger.info(msg)
                    continue
                else:
                    raise mc.CadcException(msg)
            if wl_max < w_max:
                w_max = wl_max
            if wl_min > w_min:
                w_min = wl_min
        return w_max, w_min

    def reset_energy(self, observation_type):
        result = False
        om_filter_name = json_lookup.get('filter_name')
        if (
            observation_type in ['BIAS', 'DARK']
            or (
                self._name in [
                    Inst.GMOS,
                    Inst.GMOSN,
                    Inst.GMOSS,
                ]
                and observation_type in ['BIAS', 'MASK']
            )
            or (
                om_filter_name is not None and
                ('blank' in om_filter_name or 'Blank' in om_filter_name)
            )
            or (
                self.filter_name is not None and
                ('unknown' in self.filter_name or self.filter_name == '')
            )
        ):
            logging.info(
                f'No chunk energy for {self.obs_id} obs type '
                f'{observation_type} filter name {om_filter_name}'
            )
            # 'unknown' in filter_name test obs is GN-2004B-Q-30-15-002
            # DB 23-04-19 - GN-2004B-Q-30-15-002: no energy
            # GMOS GS-2005A-Q-26-12-001.  Lots of missing metadata, including
            # release date so no energy (filter_name == '')
            # DB 04-17-21
            # BIASes and DARKs for all instruments should ignore spectral wcs
            result = True
        return result

    def reset_position(self, headers, observation_type):
        """
        Return True if there should be no spatial WCS information created at
        the chunk level.
        """
        result = False
        types = json_lookup.get('types')
        ra = self.get_ra(headers[0])
        if (
            ('AZEL_TARGET' in types and ra is None) or
            observation_type in ['BIAS', 'DARK']
        ):
            # DB - 02-04-19 - Az-El coordinate frame likely means the
            # telescope was parked or at least not tracking so spatial
            # information is irrelevant.

            # DB - 09-04-19 - AZEL_TARGET should likely be checked for all
            # datasets, and means the spatial WCS should be ignored. since
            # this generally means the telescope is not tracking and so
            # spatial WCS info isn’t relevant since the position is changing
            # with time.

            # DB 04-17-21
            # BIASes and DARKs for all instruments should ignore spatial wcs

            result = True
        # if observation_type in ['BIAS', 'DARK']:
        #     result = True
        return result

    def update_energy(self):
        raise NotImplementedError

    def update_position(self):
        # the default is to do nothing, so not a NotImplemented exception
        pass

    def update_time(self):
        # the default is to do nothing
        pass

    @staticmethod
    def get_crpix1(header):
        return 1.0

    @staticmethod
    def get_crpix2(header):
        return 1.0

    @staticmethod
    def get_sky_coord(header, ra_key, dec_key):
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
            ra_deg, dec_deg = ac.build_ra_dec_as_deg(ra_hours, dec_hours)
        return ra_deg, dec_deg


class Bhros(InstrumentType):
    def __init__(self, name):
        super(Bhros, self).__init__(name)

    def get_dec(self, header):
        # bHROS, TEXES ra/dec not in json
        return header.get('DEC')

    def get_ra(self, header):
        # bHROS, TEXES: ra/dec not in json
        return header.get('RA')

    def update_energy(self):
        """bhros-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin update_energy {self._name}')
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
        # The approximate central wavelength is json ‘central_wavelength’
        # value. Unfortunately the wavelength coverage is not straightforward.
        # See
        # http://www.gemini.edu/sciops/instruments/hros/hrosDispersion.html.
        # The CCD did not cover the entire spectrum so only a subset of the
        # entire optical spectral region was observed.

        # Use central wavelength in microns and +/- 0.2 microns as a better
        # guess-timate rather than imply that entire spectrum is present.

        if self.data_product_type == DataProductType.SPECTRUM:
            logging.debug(
                f'{self._name}: SpectralWCS spectroscopy for {self.obs_id}.'
            )
            self.fm = svofps.FilterMetadata()
            self.fm.central_wl = json_lookup.get('central_wavelength')
            self.fm.adjust_bandpass(0.2)
            self.fm.resolving_power = 150000.0
            ccd_sum = self._headers[0].get('CCDSUM')
            if ccd_sum is not None:
                temp = float(ccd_sum.split()[1])
                self.fm.resolving_power = 150000.0 / temp
            self.build_chunk_energy()
        else:
            raise mc.CadcException(
                f'{self._name}: mystery data product type '
                f'{self.data_product_type} for {self.obs_id}'
            )
        self._logger.debug(f'End update_energy {self._name}')

    def update_position(self):
        self._headers[0]['CTYPE1'] = 'RA---TAN'
        self._headers[0]['CTYPE2'] = 'DEC--TAN'
        self._headers[0]['CUNIT1'] = 'deg'
        self._headers[0]['CUNIT2'] = 'deg'
        self._headers[0]['CRVAL1'] = self.get_ra(self._headers[0])
        self._headers[0]['CRVAL2'] = self.get_dec(self._headers[0])
        self._headers[0]['CDELT1'] = RADIUS_LOOKUP[self._name]
        self._headers[0]['CDELT2'] = RADIUS_LOOKUP[self._name]
        self._headers[0]['CROTA1'] = 0.0
        self._headers[0]['NAXIS1'] = 1
        self._headers[0]['NAXIS2'] = 1
        self._headers[0]['CRPIX1'] = InstrumentType.get_crpix1(
            self._headers[0]
        )
        self._headers[0]['CRPIX2'] = InstrumentType.get_crpix2(
            self._headers[0]
        )
        self._headers[0]['CD1_1'] = self.get_cd11(self._headers[0])
        self._headers[0]['CD1_2'] = 0.0
        self._headers[0]['CD2_1'] = 0.0
        self._headers[0]['CD2_2'] = self.get_cd22(self._headers[0])
        self._headers[0]['EQUINOX'] = mc.to_float(
            self._headers[0].get('TRKEQUIN')
        )
        wcs_parser = fits2caom2.WcsParser(
            self._headers[0], self.obs_id, self._extension)
        if self.chunk is None:
            self.chunk = Chunk()
        wcs_parser.augment_position(self.chunk)
        self.chunk.position_axis_1 = 1
        self.chunk.position_axis_2 = 2
        self.chunk.position.coordsys = self._headers[0].get('TRKFRAME')
        self._logger.debug('End _update_chunk_position')


class Cirpass(InstrumentType):
    def __init__(self, name):
        super(Cirpass, self).__init__(name)

    def get_cd11(self, header, keyword='CDELT1'):
        lens_scl = header.get('LENS_SCL')
        return mc.to_float(lens_scl) / 3600.0

    def get_cd22(self, header):
        return self.get_cd11(header)

    def get_data_product_type(self, header):
        return DataProductType.SPECTRUM

    def get_dec(self, header):
        # DB - 06-03-19 - Must use FITS header info for most WCS info
        ra, dec = InstrumentType.get_sky_coord(header, 'TEL_RA', 'TEL_DEC')
        return dec

    def get_exposure(self, header):
        # exposure_time is null in the JSON
        # DB - 06-03-19 Use value of header keyword EXP_TIME for exposure time.
        return mc.to_float(header.get('EXP_TIME'))

    def get_obs_intent(self, header):
        data_label = header.get('DATALAB')
        if data_label is not None and '-CAL' in data_label:
            result = ObservationIntentType.CALIBRATION
        else:
            obs_type = header.get('OBS_TYPE')
            if obs_type in CAL_VALUES:
                result = ObservationIntentType.CALIBRATION
            else:
                result = ObservationIntentType.SCIENCE
        return result

    def get_obs_type(self, header):
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

    def get_ra(self, header):
        # DB - 06-03-19 - Must use FITS header info for most WCS info
        ra, dec = InstrumentType.get_sky_coord(header, 'TEL_RA', 'TEL_DEC')
        return ra

    def update_energy(self):
        self._logger.debug(f'Begin update_energy {self._name}')
        self.filter_name = ''
        if self.data_product_type == DataProductType.SPECTRUM:
            logging.debug(
                f'{self._name}: SpectralWCS spectral mode for {self.obs_id}.'
            )
            self.fm = svofps.FilterMetadata()
            self.fm.set_central_wl(1.0, 1.67)
            self.fm.set_bandpass(1.0, 1.67)
            self.fm.resolving_power = 3200.0
        else:
            raise mc.CadcException(
                f'{self.name}: mystery data product type '
                f'{self.data_product_type} for {self.obs_id}'
            )
        self.build_chunk_energy()
        self._logger.debug(f'End update_energy {self._name}')

    def update_position(self):
        self._headers[0]['CTYPE1'] = 'RA---TAN'
        self._headers[0]['CTYPE2'] = 'DEC--TAN'
        self._headers[0]['CUNIT1'] = 'deg'
        self._headers[0]['CUNIT2'] = 'deg'
        self._headers[0]['CRVAL1'] = self.get_ra(self._headers[0])
        self._headers[0]['CRVAL2'] = self.get_dec(self._headers[0])
        self._headers[0]['CROTA1'] = 0.0
        # So perhaps try:
        #     NAXIS1 = 33
        #     NAXIS2 = 15
        self._headers[0]['NAXIS1'] = 33
        self._headers[0]['NAXIS2'] = 15

        # TODO TODO TODO - add the get_crpix1, get_crpix2, get_cd11, get_cd22
        # methods, because that may allow for extraction of common code??????

        # LENS_SCL determines the scale/lenslet:  0.36 or 0.25 (arcseconds
        # per lens)
        #     if 0.36 then FOV is 13.0" x 4.7" (RA and Dec)
        #     if 0.25 then FOV is 9.3" x 3.5"
        lens_scl = self._headers[0].get('LENS_SCL')
        if lens_scl == '0.36':
            self._headers[0]['CDELT1'] = 13.0 / 3600.0
            self._headers[0]['CDELT2'] = 4.7 / 3600.0
        else:
            self._headers[0]['CDELT1'] = 9.3 / 3600.0
            self._headers[0]['CDELT2'] = 3.5 / 3600.0
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
        self._headers[0]['CRPIX1'] = InstrumentType.get_crpix1(
            self._headers[0]
        )
        self._headers[0]['CRPIX2'] = InstrumentType.get_crpix2(
            self._headers[0]
        )
        self._headers[0]['CD1_1'] = self.get_cd11(self._headers[0])
        self._headers[0]['CD1_2'] = 0.0
        self._headers[0]['CD2_1'] = 0.0
        self._headers[0]['CD2_2'] = self.get_cd22(self._headers[0])
        wcs_parser = fits2caom2.WcsParser(
            self._headers[0], self.obs_id, self._extension)
        if self.chunk is None:
            self.chunk = Chunk()
        wcs_parser.augment_position(self.chunk)
        self.chunk.position_axis_1 = 1
        self.chunk.position_axis_2 = 2
        self._logger.debug('End _update_chunk_position')


class F2(InstrumentType):
    def __init__(self, name):
        super(F2, self).__init__(name)

    def get_obs_type(self, header):
        # DB 03-06-21
        # check the 'types' JSON value
        result = super(F2, self).get_obs_type(header)
        types = json_lookup.get('types')
        if 'DARK' in types:
            result = 'DARK'
        return result

    def update_energy(self):
        self._logger.debug('Begin update_energy')
        # DB - 02-05-19
        # For F2 use SVO filter service for bandpass or ‘delta’ for images and
        # spectroscopy.  Treat images as for other instruments but…
        #
        # if gemini_md[‘spectroscopy’] == ‘true’:
        #  ref_wl = gemini_md[‘central_wavelength’]   or  GRWLEN header value
        #  grism = gemini_md[‘disperser’]  or  GRISM header value
        #  if gemini_md[‘mode’] == ‘LS’:  # long-slit
        #    slit_width = MASKNAME header value, e.g. ‘4pix-slit’, but need
        #                                                          only ‘4’
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
        object_value = json_lookup.get('object')
        if 'COVER CLOSED' in object_value or 'Undefined' in self.filter_name:
            # DB 30-04-19
            # Flamingos ‘Undefined’ filter:  no spectral WCS
            reset_energy = True
        else:
            filter_md = svofps.get_filter_metadata(
                self._name,
                self.filter_name,
                json_lookup.get('telescope'),
            )
            if self.data_product_type == DataProductType.IMAGE:
                self._logger.debug(
                    f'SpectralWCS: F2 imaging mode for {self.obs_id}.'
                )
                self.fm = filter_md
            elif self.data_product_type == DataProductType.SPECTRUM:
                self._logger.debug(
                    f'SpectralWCS: F2 LS|Spectroscopy mode for {self.obs_id}.'
                )
                fp_mask = self._headers[0].get('MASKNAME')
                mode = json_lookup.get('mode')
                slit_width = None
                if mode == 'LS':
                    slit_width = fp_mask[0]
                self.fm = svofps.FilterMetadata()
                self.fm.central_wl = filter_md.central_wl
                self.fm.bandpass = filter_md.bandpass
                grism_name = self._headers[0].get('GRISM')
                self._logger.debug(
                    f'F2: grism name is {grism_name} fp_mask is {fp_mask} for '
                    f'{self.obs_id}'
                )
                # lookup values from
                # https://www.gemini.edu/sciops/instruments/flamingos2/spectroscopy/longslit-spectroscopy
                lookup = {
                    '1': [1300.0, 3600.0],
                    '2': [900.0, 2800.0],
                    '3': [600.0, 1600.0],
                    '4': [350.0, 1300.0],
                    '6': [130.0, 1000.0],
                    '8': [100.0, 750.0],
                }
                if slit_width is None or slit_width not in lookup:
                    # DB 02-04-19
                    # For F2 at line 1409 of main_app.py set slit_width = ‘2’
                    # as a default of slit_width[0] is not a numeric value
                    slit_width = '2'
                if grism_name.startswith('R3K_'):
                    self.fm.resolving_power = lookup[slit_width][1]
                else:
                    self.fm.resolving_power = lookup[slit_width][0]
            else:
                raise mc.CadcException(
                    f'F2: Do not understand DataProductType '
                    f'{self.data_product_type} for {self.obs_id}'
                )

        if reset_energy:
            self._logger.info(
                f'Setting spectral WCs to none for {self.obs_id}'
            )
            cc.reset_energy(self.chunk)
        else:
            self.build_chunk_energy()
        self._logger.debug('End _update_energy')

    def update_time(self):
        """F2 FITS files have a CD3_3 element that's not supported by
        fits2caom2, so using the blueprint will not work to adjust that
        value. Set delta specifically here."""
        self._logger.debug(f'Begin update_time {self.obs_id}')
        mc.check_param(self.chunk, Chunk)
        if (
                self.chunk.time is not None
                and self.chunk.time.axis is not None
                and self.chunk.time.axis.function is not None
        ):
            exposure = mc.to_float(json_lookup.get('exposure_time'))
            self.chunk.time.axis.function.delta = mc.convert_to_days(exposure)
            logging.info(f'Updated time delta for {self.obs_id}')
        self._logger.debug(f'End update_time {self.obs_id}')


class Flamingos(InstrumentType):
    def __init__(self, name):
        super(Flamingos, self).__init__(name)

    def get_art_product_type(self, header):
        # DB - 28-02-19
        # Another relatively minor thing for FLAMINGOS:  get_obs_intent likely
        # has to have FLAMINGOS added to it.  If obs_type is DARK,
        # ACQUISIITON, FLAT or ARC then the intent is calibration, otherwise
        # science. cal_values should have those obs_type values in it as well.
        # So if ‘CAL’ is not in the data_label and obs_type is not in
        # cal_values then it’s science.
        object_value = json_lookup.get('object')
        for ii in ['Dark', 'DARK', 'Arc', 'flat', 'Flat', 'Acq', 'acq']:
            if ii in object_value:
                return ProductType.CALIBRATION
        self._get_data_label(header)
        if self._data_label is not None and 'CAL' in self._data_label:
            return ProductType.CALIBRATION
        return ProductType.SCIENCE

    def get_data_product_type(self, header):
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
        if decker is None:
            raise mc.CadcException(
                f'No mode information found for '
                f'{json_lookup.get("filename")}'
            )
        else:
            if decker == 'imaging':
                data_type = DataProductType.IMAGE
            elif decker in ['mos', 'slit']:
                grism = header.get('GRISM')
                if grism is not None:
                    if 'open' in grism:
                        data_type = DataProductType.IMAGE
                    elif 'dark' in grism:
                        data_type = DataProductType.SPECTRUM
        return data_type

    def get_exposure(self, header):
        # exposure_time is null in the JSON
        # DB - 06-03-19 Use value of header keyword EXP_TIME for exposure
        # time.
        return mc.to_float(header.get('EXP_TIME'))

    def get_obs_type(self, header):
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
                f'No mode information found for '
                f'{json_lookup.get("filename")}'
            )
        if obs_type is None:
            # DB - Also, since I’ve found FLAMINGOS spectra if OBJECT keyword
            # or json value contains ‘arc’ (any case) then it’s an ARC
            # observation type
            #
            # For FLAMINGOS since OBS_TYPE seems to always be set to ‘Object’
            # could in principal look at the OBJECT keyword value or json
            # value: if it contains ‘flat’ as a substring (any case) then set
            # observation type to ‘flat’.  Ditto for ‘dark’
            object_value = json_lookup.get('object')
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
        return obs_type

    def get_time_function_val(self, header):
        # Another FLAMINGOS correction needed:  DATE-OBS in header and
        # json doesn’t include the time  but sets it to 00:00:00.  You have
        # two choices:  concatenate DATE-OBS and UTC header values to the
        # standard form “2002-11-06T07:06:00.3” and use as you use json
        # value for computing temporal WCS data or, for FLAMINGOS only,
        # use the MJD header value in the header as the CRVAL for time.
        return header.get('MJD')

    def reset_position(self, headers, observation_type):
        result = super(Flamingos, self).reset_position(
            headers, observation_type
        )
        ra_tel = headers[0].get('RA_TEL')
        if ra_tel == 'Unavailable':
            result = True
        return result

    def update_energy(self):
        """Flamingos-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin update_energy {self._name}')

        # DB - 18-02-19 - FLAMINGOS spectral WCS should be similar to what was
        # done for NIRI.  Use NAXIS1 keyword value to determine number of
        # pixels.  The GRISM and FILTER keywords give the same filter ID.
        # The SVO filter information for KPNO/Flamingos has a blue leak in the
        # JK blocking filter which give too large a spectral range.  Can you
        # instead hard code min/max wavelengths using the top two lines in
        # this table 7:
        # http://www-kpno.kpno.noao.edu/manuals/flmn/flmn.user.html#flamspec.
        # Use these for min/max wavelengths and use the average as the
        # ‘central’ wavelength.  Spectral resolution is fixed at about 1300
        # for both grisms.

        # spectral wcs units are microns, values from Table 7 are angstroms.
        # The conversion is here.

        # DB - 21-02-19 - JH central=1.45 um, FWHM=0.95 um
        # 0 = central wavelength
        # 1 = FWHM
        lookup = {
            'JH': [1.45, 0.95],
            'HK': [(2.7588 + 1.0347) / 2.0, (2.7588 - 1.0347)],
        }
        self.fm = svofps.FilterMetadata()
        if self.filter_name in lookup:
            self.fm.central_wl = lookup[self.filter_name][0]
            self.fm.bandpass = lookup[self.filter_name][1]
        else:
            self.fm = svofps.get_filter_metadata(
                self._name,
                self.filter_name,
                json_lookup.get('telescope'),
            )
            if self.fm is None:
                raise mc.CadcException(
                    f'{self._name}: Mystery filter {self.filter_name} for '
                    f'{self.obs_id}'
                )
        if self.data_product_type == DataProductType.SPECTRUM:
            logging.debug(f'{self._name}: SpectralWCS for {self.obs_id}.')
            self.fm.resolving_power = 1300.0
        elif self.data_product_type == DataProductType.IMAGE:
            logging.debug(
                f'{self._name}: SpectralWCS imaging mode for {self.obs_id}.'
            )
        else:
            raise mc.CadcException(
                f'{self._name}: mystery data product type '
                f'{self.data_product_type} for {self.obs_id}'
            )
        self.build_chunk_energy()
        self._logger.debug(f'Begin update_energy {self._name}')

    def update_position(self):
        # DB - I see nothing in astropy that will do a transformation from
        # crota form to CD matrix, but this is it:

        # cd1_1 = cdelt1 * cos (crota1)
        # cd1_2 = -cdelt2 * sin (crota1)
        # cd2_1 = cdelt1 * sin (crota1)
        # cd2_2 = cdelt2 * cos (crota1)

        # Note that there is not a crota2 keyword (it would have the same
        # value as crota1 if it existed)
        if (
                self.chunk is not None
                and self.chunk.position is not None
                and self.chunk.position.axis is not None
                and self.chunk.position.axis.function is not None
        ):
            header = self._headers[self._extension]
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
                    self.chunk.position.axis.function.cd11 = c_delt1 * math.cos(c_rota1)
                    self.chunk.position.axis.function.cd12 = -c_delt2 * math.sin(
                        c_rota1
                    )
                    self.chunk.position.axis.function.cd21 = c_delt1 * math.sin(c_rota1)
                    self.chunk.position.axis.function.cd22 = c_delt2 * math.cos(c_rota1)
                else:
                    self._logger.info(
                        f'Missing spatial wcs inputs for {self.obs_id}'
                    )
                    self.chunk.position.axis.function.cd11 = None
                    self.chunk.position.axis.function.cd12 = None
                    self.chunk.position.axis.function.cd21 = None
                    self.chunk.position.axis.function.cd22 = None
            else:
                # DB 04-12-19
                # FLAMINGOS GS-CAL20020623-14-0080 02jun23.0080.fits
                # The header has “CRVAL1  =           3581.13808 ” which is
                # supposed to be in degrees and shouldn’t be > 360. Skip
                # spatial WCS if errors like this occur.
                self._logger.warning(
                    f'Spatial WCS set to None for {self.obs_id} '
                    f'because CRVAL1 == {crval1} and CRVAL2 == {crval2}.'
                )
                cc.reset_position(self.chunk)
        else:
            self._logger.info(f'Missing spatial wcs for {self.obs_id}')


class Fox(InstrumentType):
    def __init__(self, name):
        super(Fox, self).__init__(name)

    def get_data_product_type(self, header):
        # DB 31-08-20
        # Both cameras are used for speckle imaging.  Datasets consist of
        # cubes of 256 x 256 x 1000 images.  i.e. 1000 short exposures << 1
        # second long with 256 x 256 pixel images (or smaller images if the
        # CCD is binned).
        #
        # So dataProductType = cube
        #
        # PD 02-09-20
        # if those two files are images in the normal sense then it could make
        # sense to create separate planes with dataProductType = image that
        # end up with the correct (distinct) energy metadata.
        return DataProductType.IMAGE

    def get_target_type(self):
        return TargetType.OBJECT

    def update_energy(self):
        """General chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin update_energy {self._name}')
        if self.data_product_type is DataProductType.IMAGE:
            logging.debug(
                f'{self._name} Spectral WCS {self.data_product_type} mode for '
                f'{self.obs_id}.'
            )
            self.fm = svofps.get_filter_metadata(
                self._name,
                self.filter_name,
                json_lookup.get('telescope'),
            )
            if self.fm is None:
                raise mc.CadcException(
                    f'{self._name}: mystery filter {self.filter_name}'
                )
            self.build_chunk_energy()
        else:
            raise mc.CadcException(
                f'{self._name} no Spectral WCS support when DataProductType '
                f'{self.data_product_type} for {self.obs_id}'
            )
        self._logger.debug(f'End update_energy {self._name}')

    def update_position(self):
        self._logger.debug(f'Begin update_position for {self.obs_id}')
        if (
                self.chunk is not None
                and self.chunk.position is not None
                and self.chunk.position.axis is not None
                and self.chunk.position.axis.axis1.ctype == 'RA--TAN'
        ):
            self.chunk.position.axis.axis1.ctype = 'RA---TAN'
        self._logger.debug('End update_position.')

    def update_time(self):
        """
        DB 02-09-20
        Exposure time using JSON values isn’t correct.  I know that for this
        example Gemini shows the exposure time is 0.02 seconds but there are
        1000 x 0.02-second exposures in the cube.  The keyword EXPOSURE gives
        the total exposure time (in seconds), time.exposure, or 20 in this
        case while the json exposure_time should be the time.resolution.
        """
        self._logger.debug(f'Begin _update_chunk_time_fox {self.obs_id}')
        if self.chunk.time is not None:
            self.chunk.time.exposure = self._headers[self._extension].get(
                'EXPOSURE'
            )
            # chunk.time.resolution already set by blueprint
        self._logger.debug(f'End _update_chunk_time_fox {self.obs_id}')


class Gnirs(InstrumentType):
    def __init__(self, name):
        super(Gnirs, self).__init__(name)

    def get_data_product_type(self, header):
        # DB - 03-04-19
        # if disperser for GNIRS = MIRROR
        # then it’s an image, otherwise a spectrum.
        result = DataProductType.SPECTRUM
        if json_lookup.get('disperser') == 'MIRROR':
            result = DataProductType.IMAGE
        return result

    def reset_energy(self, observation_type):
        # DB 16-04-19
        # Energy WCS should be ignored for ‘Moving’ since we don’t know
        # what might have been in the light path when the exposure was
        # actually being taken.
        #
        # 'Dark' test obs is GN-2013B-Q-93-147-036
        # DB 13-05-19
        # GNIRS “Dark” trumps the “L” filter, so no energy.
        return (
            super(Gnirs, self).reset_energy(observation_type) or (
                self.filter_name is not None and (
                    'Moving' in self.filter_name or
                    'Dark' in self.filter_name or
                    'DARK' in self.filter_name
                )
            )
        )

    def update_energy(self):
        self._logger.debug('Begin update_energy')

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
        # would be (1.91+2.49)/2.0 or 2.2 micro^ns.  The way it is now with a
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

        # DB 10-07-19
        # Info is available here:
        # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/orderblocking-acq-nd-filters.
        # This filter is used for acquisition of targets in spectroscopy mode
        # sometimes.  Use the central wavelength of 3.295 microns and width of
        # 1.5% of that or 0.049425 microns.

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

        # DB 21-08-19
        # GN-2010B-SV-142-300-002:  add a line in the ‘10’ section of
        # gnirs_lookup similar to the ‘L’ entry but with the PAH filter
        # bandpass:  ‘PAH’: [3.270, 3.320, 570, 1800].  These are likely
        # mis-identified acquisition observations rather than science.

        # DB 31-03-21
        # It looks like Gemini changed FITS header content for gratings, etc. at
        # some point in time.  There are two prisms used in cross-dispersed mode
        # and these are identified by SXD and LXD but the leading letter.   If
        # you find an “&XD” entry in the JSON disperser value, then look at the
        # value of the FITS header PRISM = 'SXD_G5509' to find the letter that
        # should be in front of the ‘XD’ in the JSON result? It might not be the
        # first letter in the value: new files have PRISM   = 'SB+SXD_G5536'

        gnirs_lookup = {
            '10': {
                'X': [1.03, 1.17, 570, 2100],
                'J': [1.17, 1.37, 570, 1600],
                'H': [1.47, 1.80, 570, 1700],
                'K': [1.91, 2.49, 570, 1700],
                'L': [2.80, 4.20, 570, 1800],
                'M': [4.40, 6.00, 570, 1200],
                'LB+LXD': [0.9, 2.5, 1800],
                'LB+SXD': [1.2, 2.5, 1800],
                'PAH': [3.270, 3.320, 570, 1800],
            },
            '32': {
                'X': [1.03, 1.17, 1700, 5100, 1400],
                'J': [1.17, 1.37, 1600, 4800, 1400],
                'H': [1.49, 1.80, 1700, 5100, 1400],
                'K': [1.91, 2.49, 1700, 5100, 1300],
                'L': [2.80, 4.20, 1800, 5400, 1800],
                'M': [4.40, 6.00, 1240, 3700, 1240],
                'SB+SXD': [0.9, 2.5, 1800, 1800],
                'SB+LXD': [0.9, 2.5, 5400, 5400],
                'LB+LXD': [0.9, 2.5, 5400, 5400],
                'LB+SXD': [1.2, 2.5, 1800, 1800],
                'PAH': [3.270, 3.320, None, None],
            },
            '111': {
                'X': [1.03, 1.17, 6600, 17800],
                'J': [1.17, 1.37, 7200, 17000],
                'H': [1.49, 1.80, 5900, 17800],
                'K': [1.91, 2.49, 5900, 17800],
                'L': [2.80, 4.20, 6400, 19000],
                'M': [4.40, 6.00, 4300, 12800],
                'SB+SXD': [0.9, 2.5, 5400],
                'LB+SXD': [0.9, 2.5, 5400],
                'LB+LXD': [0.9, 2.5, 17000],
                'PAH': [3.270, 3.320, None, None],
            },
        }

        reset_energy = False
        if 'open' in self.filter_name:
            # DB 21-08-19
            # GS-2004A-DD-8-10-001:  set the bandpass to between 1 micron and
            # 6 microns, same as for GMOS imaging with ‘open’ filter.  This is
            # just an acquisition observation.
            w_min = 1.0
            w_max = 6.0
            self.fm = svofps.FilterMetadata()
            self.fm.set_bandpass(w_max, w_min)
            self.fm.set_central_wl(w_max, w_min)
            self.fm.set_resolving_power(w_max, w_min)
        else:
            self.fm = svofps.FilterMetadata('GNIRS')
            if self.data_product_type == DataProductType.SPECTRUM:
                self._logger.info(
                    f'SpectralWCS Spectroscopy mode for {self.obs_id}.'
                )
                disperser = json_lookup.get('disperser')
                if disperser is None:
                    self._logger.info(
                        f'No disperser. No energy for {self.obs_id}'
                    )
                    reset_energy = True
                else:
                    grating = disperser.split('_')[0]
                    # 'UNKNOWN' in grating test obs GS-CAL20040924-6-006
                    # 'ENG -' in grating test obs GN-CAL20130813-22-010
                    if (
                            'UNKNOWN' in grating
                            or 'Moving' in grating
                            or 'ENG' in grating
                    ):
                        # DB 23-04-19 - GNIRS grating UNKNOWN:  no energy
                        # DB 18-04-19
                        # If grating is moving then the observation is almost
                        # certainly garbage so ignore energy WCS.
                        # DB 13-05-19
                        # No energy for the invalid “ENG - 170000" grating
                        # entry.
                        self._logger.warning(
                            f'grating is {grating}. No energy for '
                            f'{self.obs_id}'
                        )
                        reset_energy = True
                    else:
                        if grating not in gnirs_lookup:
                            raise mc.CadcException(
                                f'Mystery grating {grating} for {self.obs_id}'
                            )

                        camera = json_lookup.get('camera')
                        if camera == 'Moving':
                            # DB 21-08-19
                            # No energy for the ‘moving’ camera
                            self._logger.info(
                                f'Camera is moving. No energy for '
                                f'{self.obs_id}'
                            )
                            reset_energy = True
                        else:
                            focal_plane_mask = json_lookup.get(
                                'focal_plane_mask'
                            )
                            if (
                                    focal_plane_mask is None
                                    or 'arcsec' not in focal_plane_mask
                            ):
                                # DB 24-04-19
                                # Assume slit width of 1.0 for GNIRS
                                # observation without a focal plane mask
                                # value.
                                slit_width = 1.0
                            else:
                                slit_width = mc.to_float(
                                    focal_plane_mask.split('arcsec')[0]
                                )

                            if 'XD' in disperser:
                                self._logger.debug(
                                    f'cross dispersed mode for {self.obs_id}.'
                                )
                                # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/crossdispersed-xd-spectroscopy/xd-prisms
                                #
                                # 0 = lower
                                # 1 = upper
                                # 2 = spectral resolution with 2-pix wide slit
                                # DB - 04-03-19 - Change the last number in
                                # each row to 1800.0 since the resolving power
                                # is the same for all 3 cases
                                #
                                # Add line to find grating ID (from long-slit
                                # code):
                                # grating = disperser.split(‘_’)[0]
                                #
                                # Then change ‘lookup’ to include grating.
                                #
                                # I can’t find any other combinations (e.g.
                                # ‘LB+LXD+32’) but no guarantee that I won’t
                                # have to add another line or two if we see
                                # failures. Wavelength coverage isn’t correct
                                # for the R=5400 entries because only about
                                # 1/3rd of the full band pass is covered but
                                # in bits and pieces that we can’t identify in
                                # raw image.

                                lookup = None
                                if '&XD' in disperser:
                                    prism = self._headers[0].get('PRISM')
                                    cd_mode = prism[prism.index('XD') - 1]
                                    coverage = f'{cd_mode}XD'
                                else:
                                    coverage = disperser[-3:]

                                if camera.startswith('Short'):
                                    lookup = f'SB+{coverage}'
                                elif camera.startswith('Long'):
                                    lookup = f'LB+{coverage}'

                                if camera.startswith('Long'):
                                    slit_table_value = 0.1
                                    lookup_index = 2
                                elif camera.startswith('Short'):
                                    slit_table_value = 0.3
                                    lookup_index = 2
                                else:
                                    raise mc.CadcException(
                                        f'Mystery camera definition {camera} '
                                        f'for {self.obs_id}'
                                    )
                            else:
                                self._logger.debug(
                                    f'long slit mode for {self.obs_id}.'
                                )
                                if self.filter_name == 'PAH':
                                    lookup = self.filter_name
                                else:
                                    bandpass = self.filter_name[0]
                                    lookup = bandpass
                                    if camera.startswith('Long'):
                                        slit_table_value = 0.1
                                        lookup_index = 3
                                    elif camera.startswith('Short'):
                                        date_time = ac.get_datetime(
                                            json_lookup.get(
                                                'ut_datetime'
                                            )
                                        )
                                        if date_time > ac.get_datetime(
                                                '2012-11-01T00:00:00'
                                        ):
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
                                            f'GNIRS: Mystery camera '
                                            f'definition {camera} for '
                                            f'{self.obs_id}'
                                        )

                            lookup = lookup.upper()
                            if lookup not in gnirs_lookup[grating]:
                                raise mc.CadcException(
                                    f'GNIRS: Mystery lookup {lookup} for '
                                    f'grating {grating}, obs {self.obs_id}'
                                )
                            bounds = gnirs_lookup[grating][lookup]
                            self.fm.set_bandpass(bounds[1], bounds[0])
                            if lookup == 'PAH':
                                self.fm.resolving_power = None
                            else:
                                self.fm.resolving_power = (
                                        slit_table_value
                                        * bounds[lookup_index]
                                        / slit_width
                                )
                            self.fm.set_central_wl(bounds[1], bounds[0])
            elif self.data_product_type == DataProductType.IMAGE:
                self._logger.debug(
                    f'SpectralWCS imaging mode for {self.obs_id}.'
                )
                # https://www.gemini.edu/sciops/instruments/gnirs/imaging
                # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/orderblocking-acq-nd-filters

                # DB 23-04-19
                # GNIRS filter “L”: add lower/upper bandpass info for the
                # L and M filter on this page,
                # https://www.gemini.edu/sciops/instruments/gnirs/spectroscopy/orderblocking-acq-nd-filters
                # to ‘imaging’ dictionary.  (In ‘imaging’ I think the K order
                # blocking filter lower bandpass should be 1.91).
                # Despite Gemini’s footnote stating that L and M filters are
                # NOT used for acquisition images.

                # 0 - min wavelength
                # 1 - max wavelength
                imaging = {
                    'Y': [0.97, 1.07],
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
                    'L': [2.8, 4.2],
                }

                # DB 30-04-19
                # GN-CAL20180607-27-001:  ignore ‘PupilViewer’ if possible.
                # DB 21-08-19
                # PV+K: ‘PV’ = pupil viewer.  Shouldn’t affect transmission
                # values (I think) so ignore when looking up lower/upper
                # bandpasses.
                self.filter_name = self.filter_name.replace(
                    'PupilViewer', ''
                ).replace('PV', '')
                self.filter_name = self.filter_name.strip('+')

                if (
                        len(self.filter_name) == 1
                        or self.filter_name == 'H2'
                        or self.filter_name == 'PAH'
                        or self.filter_name == 'XD'
                ):
                    bandpass = self.filter_name
                else:
                    bandpass = self.filter_name[0]

                bandpass = bandpass.upper()
                if bandpass in imaging:
                    bounds = imaging[bandpass]
                else:
                    raise mc.CadcException(
                        f'GNIRS: Unexpected filter_name {self.filter_name} '
                        f'for {self.obs_id}'
                    )

                self.fm.set_central_wl(bounds[1], bounds[0])
                self.fm.set_bandpass(bounds[1], bounds[0])
            else:
                raise mc.CadcException(
                    f'GNIRS: Unexpected DataProductType '
                    f'{self.data_product_type} for {self.obs_id}'
                )

        if reset_energy:
            cc.reset_energy(self.chunk)
        else:
            self.build_chunk_energy()
        self._logger.debug('End update_energy')


class Gpi(InstrumentType):
    def __init__(self, name):
        super(Gpi, self).__init__(name)

    def get_cd11(self, header):
        return RADIUS_LOOKUP[self._name] / header.get('NAXIS1')

    def get_cd22(self, header):
        return self.get_cd11(header)

    def get_data_product_type(self, header):
        mode = json_lookup.get('mode')
        if mode is None:
            raise mc.CadcException(
                f'No mode information found for '
                f'{json_lookup.get("filename")}'
            )
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
                f'Mystery GPI mode {mode} for '
                f'{json_lookup.get("filename")}'
            )
        return result

    def update_energy(self):
        self._logger.debug('Begin _update_chunk_energy_gpi')

        # DB - 22-02-19
        # For imaging energy WCS, use standard imaging algorithm. i.e central
        # wavelength, bandpass from SVO filter, and
        # resolution = central_wavelength/bandpass
        #
        # For energy WCS:
        #
        # naxis = extension header NAXIS2 value
        #
        # Need to hard-code the resolution from values in top table on this
        # page:
        # https://www.gemini.edu/sciops/instruments/gpi/observing-modes-metamodes.
        # Use mean of column 4 range.  e.g. for H filter use 46.5
        # units are microns
        gpi_lookup = {
            'Y': (36 + 34) / 2.0,
            'J': (39 + 35) / 2.0,
            'H': (49 + 44) / 2.0,
            'K1': (70 + 62) / 2.0,
            'K2': (83 + 75) / 2.0,
        }

        filter_md = svofps.get_filter_metadata(
            self._name,
            self.filter_name,
            json_lookup.get('telescope'),
        )
        if self.data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'SpectralWCS imaging mode for {self.obs_id}.')
            self.fm = filter_md
        elif self.data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(
                f'SpectralWCS Spectroscopy mode for {self.obs_id}.'
            )
            self.fm = svofps.FilterMetadata()
            self.fm.central_wl = filter_md.central_wl
            self.fm.bandpass = filter_md.bandpass
            if self.filter_name in gpi_lookup:
                self.fm.resolving_power = gpi_lookup[self.filter_name]
            else:
                raise mc.CadcException(
                    f'GPI: Mystery filter name {self.filter_name} for '
                    f'resolving power {self.obs_id}'
                )
        else:
            raise mc.CadcException(
                f'GPI: Do not understand DataProductType '
                f'{self.data_product_type} for {self.obs_id}'
            )

        self.build_chunk_energy()
        self._logger.debug('End update_energy')

    def update_position(self):
        self._logger.debug('Begin update_position')

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

        header = self._headers[1]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(header)
        header['CRVAL2'] = self.get_dec(header)
        header['CRPIX1'] = Gpi.get_crpix1(header)
        header['CRPIX2'] = Gpi.get_crpix2(header)
        header['CD1_1'] = self.get_cd11(header)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(header)

        wcs_parser = fits2caom2.WcsParser(
            header, self.obs_id, self._extension
        )
        if self.chunk is None:
            self.chunk = Chunk()
        wcs_parser.augment_position(self.chunk)
        self.chunk.position_axis_1 = 1
        self.chunk.position_axis_2 = 2
        self.chunk.position.coordsys = header.get('RADESYS')
        if self._extension == 1:
            # equinox information only available from
            # 0th header
            equinox = self._headers[0].get('TRKEQUIN')
            if (
                equinox is not None
                and 1800.0 <= equinox <= 2500.0
            ):
                self.chunk.position.equinox = equinox
            else:
                # DB 07-06-21
                # No spatial WCS in these cases.
                cc.reset_position(self.chunk)
        self._logger.debug('End update_position')

    @staticmethod
    def get_crpix1(header, keyword='NAXIS1'):
        naxis1 = header.get(keyword)
        if naxis1 is None:
            result = None
        else:
            result = naxis1 / 2.0
        return result

    @staticmethod
    def get_crpix2(header):
        return Gpi.get_crpix1(header, 'NAXIS2')


class Graces(InstrumentType):
    def __init__(self, name):
        super(Graces, self).__init__(name)

    def get_filter_name(self):
        self._filter_name = None

    def reset_position(self, headers, observation_type):
        result = super(Graces, self).reset_position(
            headers, observation_type
        )
        # DB 23-04-19
        # Ignore spatial WCS for the GRACES dataset with EPOCH=0.0.  Not
        # important for a bias. For GMOS we skip spatial WCS for biases
        # (and maybe for some other instruments).

        # DB 24-04-19
        # GRACES:  you can ignore spatial WCS for flats if RA/Dec are not
        # available.   Ditto for GNIRS darks.

        # DB 30-04-19
        # Ignore spatial WCS for any GRACES arcs without RA/Dec values.

        ra = json_lookup.get('ra')
        dec = json_lookup.get('dec')
        if (
            (ra is None and dec is None) or
            observation_type in ['BIAS', 'FLAT', 'ARC']
        ):
            result = True
        return result

    def update_energy(self):
        """graces-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin update_energy {self._name}')

        # DB - 22-02-19  Axis 2 is the spectral axis.   Resolution is
        # 67,500 divided by second number in json ‘detector_binning’
        # value of ‘N x N’.  The central wavelength in json
        # is always set to 0.7 microns and range is from about 0.4 to
        # 1.0 microns.  So use 0.7 as the crval, crpix = naxis2/2.0 and
        # delta = (1.0 - 0.4)/naxis2.  Just a kludge to get spectral
        # range correct for an echelle spectrograph.

        if self.data_product_type == DataProductType.SPECTRUM:
            logging.debug(
                f'{self._name}: SpectralWCS {self.data_product_type} for '
                f'{self.obs_id}.'
            )
            self.fm = svofps.FilterMetadata()
            self.fm.central_wl = json_lookup.get('central_wavelength')
            self.fm.set_bandpass(1.0, 0.4)
            ccd_sum = json_lookup.get('detector_binning')
            self.fm.resolving_power = 67500.0
            if ccd_sum is not None:
                temp = float(ccd_sum.split('x')[1])
                self.fm.resolving_power = 67500.0 / temp
            self.build_chunk_energy()
        else:
            raise mc.CadcException(
                f'{self._name}: mystery data product type '
                f'{self.data_product_type} for {self.obs_id}'
            )
        self._logger.debug(f'End update_energy {self._name}')

    def update_position(self):
        if self.mode is not None and self.mode != 'imaging':
            header = self._headers[self._extension]
            header['CTYPE1'] = 'RA---TAN'
            header['CTYPE2'] = 'DEC--TAN'
            header['CUNIT1'] = 'deg'
            header['CUNIT2'] = 'deg'
            header['CRVAL1'] = self.get_ra(header)
            header['CRVAL2'] = self.get_dec(header)
            header['CDELT1'] = RADIUS_LOOKUP[self._name]
            header['CDELT2'] = RADIUS_LOOKUP[self._name]
            header['CROTA1'] = 0.0
            header['NAXIS1'] = 1
            header['NAXIS2'] = 1
            header['CRPIX1'] = InstrumentType.get_crpix1(header)
            header['CRPIX2'] = InstrumentType.get_crpix2(header)
            header['CD1_1'] = self.get_cd11(header)
            header['CD1_2'] = 0.0
            header['CD2_1'] = 0.0
            header['CD2_2'] = self.get_cd22(header)
            wcs_parser = fits2caom2.WcsParser(
                header, self.obs_id, self._extension
            )
            if self.chunk is None:
                self.chunk = Chunk()
            wcs_parser.augment_position(self.chunk)
            self.chunk.position_axis_1 = 1
            self.chunk.position_axis_2 = 2
        self._logger.debug('End update_position')


class Gsaoi(InstrumentType):
    def __init__(self, name):
        super(Gsaoi, self).__init__(name)

    def reset_energy(self, observation_type):
        # DB 24-04-19
        # ‘Unknown+Blocked2’ filter, no spectral WCS.
        return (
            super(Gsaoi, self).reset_energy(observation_type) or (
                self.filter_name is not None and (
                    'Unknown' in self.filter_name or
                    'Blocked' in self.filter_name
                )
            )
        )

    def update_energy(self):
        """General chunk-level Energy WCS construction.
        Same as Fox.
        """
        self._logger.debug(f'Begin update_energy {self._name}')
        if self.data_product_type is DataProductType.IMAGE:
            logging.debug(
                f'{self._name} Spectral WCS {self.data_product_type} mode for '
                f'{self.obs_id}.'
            )
            self.fm = svofps.get_filter_metadata(
                self._name,
                self.filter_name,
                json_lookup.get('telescope'),
            )
            if self.fm is None:
                raise mc.CadcException(
                    f'{self._name}: mystery filter {self.filter_name}'
                )
            self.build_chunk_energy()
        else:
            raise mc.CadcException(
                f'{self._name} no Spectral WCS support when DataProductType '
                f'{self.data_product_type} for {self.obs_id}'
            )
        self._logger.debug(f'End update_energy {self._name}')


class Hokupaa(InstrumentType):
    def __init__(self, name):
        super(Hokupaa, self).__init__(name)

    def get_art_product_type(self, header):
        result = ProductType.SCIENCE
        image_type = header.get('IMAGETYP')
        # took the list from the AS GEMINI Obs. Type values
        if image_type in ['DARK', 'BIAS', 'CAL', 'ARC', 'FLAT']:
            result = ProductType.CALIBRATION
        return result

    def get_cd11(self, header):
        pix_scale = header.get('PIXSCALE')
        if pix_scale is None:
            result = None
        else:
            result = mc.to_float(pix_scale) / 3600.0
        return result

    def get_cd22(self, header):
        return self.get_cd11(header)

    def get_dec(self, header):
        ra_ignore, dec = InstrumentType.get_sky_coord(header, 'RA', 'DEC')
        return dec

    def get_obs_type(self, header):
        # DB - 01-18-19 Use IMAGETYP as the keyword
        return header.get('IMAGETYP')

    def get_ra(self, header):
        ra, dec_ignore = InstrumentType.get_sky_coord(header, 'RA', 'DEC')
        return ra

    def reset_energy(self, observation_type):
        # DB 14-08-19
        # no energy since ‘home’ might either be an open position or a blocked
        # position.
        return (
            super(Hokupaa, self).reset_energy(observation_type) or (
                self.filter_name is not None and (
                    'LowFlx' in self.filter_name or 'home' in self.filter_name
                )
            )
        )

    def reset_position(self, headers, observation_type):
        result = super(Hokupaa, self).reset_position(
            headers, observation_type
        )
        ra = self.get_ra(headers[0])
        if ra is None:
            result = True
        return result

    def update_energy(self):
        """hokupaa-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin update_energy {self._name}')

        # DB - 18-01-19 - Note: it is always imaging for Hokupa'a + QUIRC so
        # that could be hardcoded.

        # DB - 20-02-19 - This page has a table of QUIRC (the camera part of
        # Hokupa’a + QUIRC) filter names and bandpasses (in microns).
        # Use this as a lookup for central wavelengths and bandpass
        # http://www.gemini.edu/sciops/instruments/uhaos/uhaosQuirc.html

        # J+CO is "Dark Position", units are microns
        # 0 - central wavelength
        # 1 - bandpass
        hokupaa_lookup = {
            'J': [1.25, 0.171],
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
            'CO': [2.29, 0.0020],
        }

        # DB 27-05-19
        # The bandpasses for the 7 bottom entries have to be corrected
        # because of the misleading column heading.  The bandpasses for these
        # 7 are in Angstroms, not microns.  So divide the values by 10,000 to
        # convert to microns. 170 Å = 0.017 microns for the FeII filter. H2/23
        # = 1-0 S(1) H2 - bandpass should be 0.0023
        # HKnotch = H+K notch
        # 1.56/120 = methane low  - bandpass should be 0.0120
        # 1.71/120 = methane high - ditto
        # LowFlx is likely the J+CO combo in line 14 of the filter table in
        # the web page:  effectively a shutter so no energy in this case. And I
        # think FeI/17 is likely meant to be the FeII filter (just going by the
        # “17” in the bandpass.  I’m not aware of an IR Fe I filter…

        # DB 21-08-19
        # GN-2001A-C-23-38-613:  filter BrG/20 -> BrG = not sure what the
        # ‘/20’ denotes since it is not the correct bandpass if that’s what
        # it’s for, but likely Br(gamma) filter, #11 on the list on the web
        # page
        #
        # GN-2001B-C-3-75-110:   filter HeI/30 -> He I], #9 on the list.  In
        # this case ‘/30’ might be the bandpass  RA/Dec values in the header
        # look reasonable but archive suggests ‘invalid’ for some reason so no
        # spatial WCS I guess.  DEC value might not agree with the telescope
        # elevation value…
        #
        # GN-2001B-C-16-19-039:  2.26/60 - it’s K-continuum.
        #
        # GN-2001A-DD-2-1-537:  Much of the header info seems to be incorrect
        # so no WCS.  But what’s odd is that searching Gemini’s archive for
        # this data label returns 3 files with different file names:  23, 24
        # and 25 February.  EQUINOX only seems to be incorrect in
        # 01FEB23_537.fits.

        filter_name_repair = {
            'H2/23': '1-0 S(1) H2',
            'HKnotch': 'H+K notch',
            '1.56/120': 'methane low',
            '1.71/120': 'methane high',
            'FeI/17': 'FeII',
            '2.26/60': 'K-continuum',
            'BrG/20': 'H Br(gamma)',
            'HeI/30': 'HeI',
        }

        if self.filter_name in filter_name_repair:
            self.filter_name = filter_name_repair[self.filter_name]
        if self.filter_name not in hokupaa_lookup:
            raise mc.CadcException(
                f'{self._name}: Mystery filter {self.filter_name} for '
                f'{self.obs_id}'
            )
        if self.data_product_type == DataProductType.IMAGE:
            logging.debug(
                f'{self._name}: SpectralWCS imaging mode for '
                f'{self.obs_id}.'
            )
            self.fm = svofps.FilterMetadata()
            self.fm.central_wl = hokupaa_lookup[self.filter_name][0]
            self.fm.bandpass = hokupaa_lookup[self.filter_name][1]
        else:
            raise mc.CadcException(
                f'{self._name}: mystery data product type '
                f'{self.data_product_type} for {self.obs_id}'
            )
        self.build_chunk_energy()
        self._logger.debug(f'End update_energy {self._name}')

    def update_position(self):
        header = self._headers[self._extension]
        equinox = mc.to_float(header.get('EQUINOX'))
        if not 1800.0 <= equinox <= 2500.0:
            self._logger.warning(
                f'EQUINOX value is wrong ({equinox}), no spatial WCS for '
                f'{self.obs_id}'
            )
            return

        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(header)
        header['CRVAL2'] = self.get_dec(header)
        header['CDELT1'] = RADIUS_LOOKUP[self._name]
        header['CDELT2'] = RADIUS_LOOKUP[self._name]
        header['CROTA1'] = 0.0
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
        header['CRPIX1'] = Gpi.get_crpix1(header, 'CRPIX1')
        header['CRPIX2'] = Gpi.get_crpix1(header, 'CRPIX2')
        header['CD1_1'] = self.get_cd11(header)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(header)

        wcs_parser = fits2caom2.WcsParser(
            header, self.obs_id, self._extension
        )
        if self.chunk is None:
            self.chunk = Chunk()
        wcs_parser.augment_position(self.chunk)
        self.chunk.position_axis_1 = 1
        self.chunk.position_axis_2 = 2
        self._logger.debug('End update_position')


class Hrwfs(InstrumentType):
    def __init__(self, name):
        super(Hrwfs, self).__init__(name)

    def update_energy(self):
        self._logger.debug('Begin _update_chunk_energy_hrwfs')
        # DB 27-05-19
        # e.g. GS-CAL20020322-7-0003 2002mar22_0055, filter1=‘neutral’,
        # filter2=‘open’, so treat similarly to “open” for GMOS since it’s a
        # similar CCD:  maybe 0.35 - 1.0 microns.
        if self.data_product_type == DataProductType.IMAGE:
            self._logger.debug(
                f'Spectral WCS {self.data_product_type} mode for '
                f'{self.obs_id}.'
            )
            if (
                'open' in self.filter_name
                or 'neutral' in self.filter_name
                or 'undefined' in self.filter_name
            ):
                w_min = 0.35
                w_max = 1.0
                self.fm = svofps.FilterMetadata()
                self.fm.set_bandpass(w_max, w_min)
                self.fm.set_central_wl(w_max, w_min)
                self.fm.set_resolving_power(w_max, w_min)
                # DB 27-05-19
                # bandpassName likely best to set to NULL
                self.filter_name = ''
            else:
                self.fm = svofps.get_filter_metadata(
                    self._name,
                    self.filter_name,
                    json_lookup.get('telescope'),
                )
                temp = []
                for ii in self.filter_name.split('+'):
                    temp.append(ii.split('_')[0])
                self.filter_name = '+'.join(i for i in temp)
            self.build_chunk_energy()
        else:
            raise mc.CadcException(
                f'hrwfs no Spectral WCS support when DataProductType '
                f'{self.data_product_type} for {self.obs_id}'
            )
        self._logger.debug('End _update_chunk_energy_hrwfs')


class Igrins(InstrumentType):
    """
    DB 21-07-21
    IGRINS metadata:
    - modify DATALAB values to give (hopefully) unique observation IDs my
    - use an observation ID of GS-2020B-Q-315-23-1104 instead of
    GS-2020B-Q-315-23-0 for file SDCH_20201104_0023.fits, by grabbing the
    MMDD from the file name (1104 in this case) and replacing the trailing
    -0 with -MMDD.

    Spatial WCS:
    - use JSON RA/Dec values for crval1 and crval2.  Same for both artifacts
    - hardcode a 5" x 5" field of view.   The spectrograph slit is
      0.34" x 5" but rotation on sky can’t be determined so best to use
      5" x 5"  square.  Same as we did for DAO 1.2m data.
    - CD11 = CD22 = 5.0/3600.0 and CD12=CD21=0.0
    - single ‘pixel’, so naxis1=naxis2=1
    - crpix1 = crpix2 = 1.0

    Temporal WCS:
    - standard things from JSON should be fine

    Energy WCS:
    - resolving power is fixed at 45,000 for both artifacts
    - JSON filter_name or wavelength_band will tell you if it’s an H- or
      K-band artifact, or the H/K in the file name.
    - H-band artifact
      - wavelength is between 1.47 and 1.83 microns
    - K-band artifact
      - wavelength is between 1.93 and 2.5 microns
   """
    def __init__(self, name):
        super(Igrins, self).__init__(name)

    def get_target_type(self):
        return TargetType.OBJECT

    def update_energy(self):
        self._logger.debug('Begin update_energy')
        # 0 minimum wavelength in microns
        # 1 max wavelength in microns
        lookup = {
            'H': [1.47, 1.83],
            'K': [1.93, 2.5]
        }
        filter_name = json_lookup.get('filter_name')
        if filter_name in lookup.keys():
            self.fm = svofps.FilterMetadata()
            self.fm.set_bandpass(
                lookup.get(filter_name)[1], lookup.get(filter_name)[0]
            )
            self.fm.set_central_wl(
                lookup.get(filter_name)[1], lookup.get(filter_name)[0]
            )
            self.fm.resolving_power = 45000.0
        else:
            raise mc.CadcException(
                f'Mystery filter {filter_name} for {self.obs_id}'
            )
        self.build_chunk_energy()
        self._logger.debug('End update_energy')

    def update_position(self):
        self._logger.debug('Begin update_position')
        header = self._headers[0]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(header)
        header['CRVAL2'] = self.get_dec(header)
        header['CDELT1'] = RADIUS_LOOKUP[self._name]
        header['CDELT2'] = RADIUS_LOOKUP[self._name]
        header['CROTA1'] = 0.0
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
        header['CRPIX1'] = 1.0
        header['CRPIX2'] = 1.0
        header['CD1_1'] = RADIUS_LOOKUP[self._name]
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = RADIUS_LOOKUP[self._name]

        wcs_parser = fits2caom2.WcsParser(
            header, self.obs_id, self._extension
        )
        if self.chunk is None:
            self.chunk = Chunk()
        wcs_parser.augment_position(self.chunk)
        self.chunk.position_axis_1 = 1
        self.chunk.position_axis_2 = 2
        self._logger.debug('End update_position')


class Michelle(InstrumentType):
    def __init__(self, name):
        super(Michelle, self).__init__(name)

    def get_filter_name(self):
        # DB - 08-04-19
        # Use header FILTERA and FILTERB values to determine the filter
        # bandpass
        # 08-04-19
        # Note that FILTERB sometimes has a value “Clear_B” and that should
        # be ignored.  Not sure if “Clear_A” is possible but maybe best to
        # allow for it.
        result = InstrumentType._search_through_keys(
            self._headers[0], ['FILTER'], FILTER_VALUES_TO_IGNORE,
        )
        self._filter_name = result

    def make_axes_consistent(self):
        super(Michelle, self).make_axes_consistent()
        if self.chunk.naxis == 4:
            if self.chunk.position is None:
                self.chunk.naxis = None
                if self.chunk.time_axis is not None:
                    self.chunk.time_axis = None
            else:
                self.chunk.naxis = 2
                self.chunk.time_axis = None

    def reset_energy(self, observation_type):
        # 'No Value' in filter_name test obs GN-2005A-C-14-45-002
        # DB 09-04-19 - “Blank-B” -> no energy.  It is a bias exposure
        # apparently and those shouldn't have energy WCS.
        # DB 23-04-19 - “No Value” -> no energy
        return (
            super(Michelle, self).reset_energy(observation_type) or
            (self.filter_name is not None and (
                'Blank' in self.filter_name or 'No Value' in self.filter_name
                )
             )
        )

    def update_energy(self):
        self._logger.debug('Begin _update_chunk_energy_michelle')

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
        lookup = {
            'LowN': [200, 2.0],
            'LowQ': [110, 3.0],
            'MedN1': [1000, 2.0],
            'MedN2': [3000, 2.0],
            'Echelle': [30000, 2.0],
        }

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
            'QBlock': [(16.1 + 25) / 2, 16.1, 25.0],
            'NBlock': [(14 + 6.8) / 2, 6.8, 14.0],
        }

        self.fm = svofps.get_filter_metadata(
            self._name,
            self.filter_name,
            json_lookup.get('telescope'),
        )
        if self.fm is None:  # means filter_name not found
            # DB 09-04-19 - Use 100 microns for the initial max for michelle.
            w_max, w_min = self.multiple_filter_lookup(michelle, wl_max=100)
            self.fm = svofps.FilterMetadata()
            self.fm.set_bandpass(w_max, w_min)
            self.fm.set_central_wl(w_max, w_min)
        if self.data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'Spectral WCS spectrum for {self.obs_id}.')
            if self.data_product_type == DataProductType.SPECTRUM:
                disperser = json_lookup.get('disperser')
                focal_plane_mask = json_lookup.get('focal_plane_mask')
                slit_width = float(focal_plane_mask.split('_')[0])
                if disperser not in lookup:
                    raise mc.CadcException(
                        f'michelle: Mystery disperser {disperser} for '
                        f'{self.obs_id}'
                    )
                self.fm.resolving_power = (
                    lookup[disperser][0] * lookup[disperser][1] / slit_width
                )
        elif self.data_product_type == DataProductType.IMAGE:
            self._logger.debug(
                f'Spectral WCS imaging mode for {self.obs_id}.'
            )
        else:
            raise mc.CadcException(
                f'michelle: no Spectral WCS support when DataProductType '
                f'{self.data_product_type} for {self.obs_id}'
            )

        # use the json value for bandpass_name value - it's representative of
        # multiple filters
        self.filter_name = json_lookup.get('filter_name')
        self.build_chunk_energy()
        self._logger.debug('End _update_chunk_energy_michelle')


class Nici(InstrumentType):
    def __init__(self, name):
        super(Nici, self).__init__(name)

    def get_filter_name(self):
        # NICI - use filter names from headers, because there's a different
        # filter/header, and the JSON summary value obfuscates that
        result = InstrumentType._search_through_keys(
            self._headers[self._extension],
            ['FILTER'],
            FILTER_VALUES_TO_IGNORE,
        )
        self._filter_name = result

    def reset_energy(self, observation_type):
        # DB 04-04-19
        # If one of the NICI filters is ‘Block’ then energy WCS should be
        # ignored for that extension.
        return (
            super(Nici, self).reset_energy(observation_type) or
            (self.filter_name is not None and 'Block' in self.filter_name)
        )

    def update_energy(self):
        self._logger.debug('Begin _update_chunk_energy_nici')
        # DB 22-08-19
        # Duplicate the line “‘FeII’: [1.644000, 1.631670, 1.656330]” in the
        # NICI filters replacing “FeII” with “[FeII]“?   Although it would be
        # better to have only [FeII] show up in the pick list and not both
        # Fe II and [Fe II]…

        # select filter_id, wavelength_central, wavelength_lower,
        #        wavelength_upper
        # from gsa..gsa_filters where instrument = 'NICI'
        # 0 - central
        # 1 - lower
        # 2 - upper
        #
        # dict with the barcodes stripped from the names as returned by query
        # DB - 22-02-19 - units are microns
        NICI = {
            'Br-gamma': [2.168600, 2.153900, 2.183300],
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
            'Mprime': [4.680000, 4.550000, 4.790000],
        }
        self.filter_name = self.filter_name.split('_G')[0]
        self.filter_name = self.filter_name.replace('[FeII]', 'FeII')
        self.fm = svofps.get_filter_metadata(
            self._name,
            self.filter_name,
            json_lookup.get('telescope'),
        )
        if self.data_product_type == DataProductType.IMAGE:
            self._logger.debug(
                f'SpectralWCS imaging mode for {self.obs_id}.'
            )
            if self.fm is None:  # means filter_name not found
                w_max, w_min = self.multiple_filter_lookup(NICI)
                self.fm = svofps.FilterMetadata()
                self.fm.set_bandpass(w_max, w_min)
                self.fm.set_central_wl(w_max, w_min)
                self.fm.set_resolving_power(w_max, w_min)

            temp = json_lookup.get('filter_name')
            # NICI has two different bandpass names (most of the time) in
            # its two chunks.  Pat says in this case nothing will be put in
            # the bandpass name for the plane.  Add code to combine the two
            # chunk bandpass names to create a plane bandpass name only for
            # this instrument
            self.filter_name = temp.replace('[FeII]', 'FeII')
            self.build_chunk_energy()
        else:
            raise mc.CadcException(
                f'NICI: Do not understand DataProductType '
                f'{self.data_product_type} from {self.obs_id}'
            )
        self._logger.debug('End _update_chunk_energy_nici')


class Nifs(InstrumentType):
    def __init__(self, name):
        super(Nifs, self).__init__(name)

    def get_cd11(self, header):
        # DB - 05-03-19 - NIFS needs a division by NAXIS1/2 for the
        # cdelta1/2 calculations.
        return RADIUS_LOOKUP[self._name] / header.get('NAXIS1')

    def get_cd22(self, header):
        return self.get_cd11(header)

    def reset_energy(self, observation_type):
        # filter_name == '', test obs is GN-CAL20050301-17-001
        # DB 18-04-19
        #
        # Gemini archive shows WaveBand=INVALID - ignore energy WCS
        # DB 23-04-19 - ‘no value’ -> no energy

        # DB 23-04-19 - Looks like NIRI observation GN-CAL20020623-1-011 is
        # skipping the ‘invalid’ values in FILTER1/2 headers.  So energy
        # should likely be skipped.  In this case the json value is
        # INVALID&INVALID.
        return (
            super(Nifs, self).reset_energy(observation_type) or
            (
                self.filter_name is not None and (
                    'Blocked' in self.filter_name or
                    'INVALID' in self.filter_name
                )
            )
        )

    def reset_position(self, headers, observation_type):
        result = super(Nifs, self).reset_position(headers, observation_type)

        # DB - 08-04-19 - json ra/dec values are null for
        # the file with things set to -9999.  Ignore
        # spatial WCS for these cases.

        # get the values from JSON directly, because the
        # function uses header values, which are set to
        # unlikely defaults
        ra = json_lookup.get('ra')
        dec = json_lookup.get('dec')
        if ra is None and dec is None:
            result = True
        elif (
                ra is not None
                and math.isclose(ra, 0.0)
                and dec is not None
                and math.isclose(dec, 0.0)
        ):
            result = True
        return result

    def update_energy(self):
        self._logger.debug('Begin _update_chunk_energy_nifs')

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

        nifs_lookup = {
            'Z': {
                'ZJ': [1.05, 0.94, 1.15, 4990.0, 60.1],
            },
            'J': {
                'ZJ': [1.25, 1.15, 1.33, 6040.0, 49.6],
                'JH': [1.25, 1.15, 1.33, 6040.0, 49.6],
            },
            'H': {
                'JH': [1.65, 1.49, 1.80, 5290.0, 56.8],
                'HK': [1.65, 1.49, 1.80, 5290.0, 56.8],
            },
            'K': {
                'HK': [2.20, 1.99, 2.40, 5290.0, 56.7],
            },
            'K_Short': {
                'HK': [2.20, 1.98, 2.41, 5290.0, 56.7],
            },
            'K_Long': {
                'HK': [2.20, 1.98, 2.41, 5290.0, 56.7],
            },
        }

        self.fm = svofps.get_filter_metadata(
            self._name,
            self.filter_name,
            json_lookup.get('telescope'),
        )

        if self.data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'spectroscopy for {self.obs_id}.')
            grating = json_lookup.get('disperser')
            if grating in nifs_lookup:
                if self.filter_name in nifs_lookup[grating]:
                    self.fm = svofps.FilterMetadata('NIFS')
                    self.fm.set_bandpass(
                        nifs_lookup[grating][self.filter_name][2],
                        nifs_lookup[grating][self.filter_name][1],
                    )
                    self.fm.central_wl = json_lookup.get(
                        'central_wavelength'
                    )
                    self.fm.resolving_power = nifs_lookup[grating][self.filter_name][3]
                    if self.fm.central_wl is None:
                        # DB 14-04-21
                        #  if no central wavelength in json then use filter
                        #  bounds
                        self.fm.central_wl = nifs_lookup[grating][self.filter_name][0]
                        self._logger.warning(
                            f'JSON central_wavelength is None for '
                            f'{self.obs_id}. Using {self.fm.central_wl} '
                            f'instead.'
                        )
                else:
                    self._logger.info(
                        f'No energy. filter_name {self.filter_name} with '
                        f'disperser {grating} for {self.obs_id}'
                    )
                    fm = None
            else:
                self._logger.info(
                    f'No energy. grating {grating} for {self.obs_id}'
                )
                fm = None
        elif self.data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'NIFS: imaging for {self.obs_id}.')
            # DB - 01-03-19
            # NIFS images should just use the standard imaging procedure
            # for resolution (central_wavelength/bandpass).
        else:
            raise mc.CadcException(
                f'NIFS: DataProductType {self.data_product_type} for '
                f'{self.obs_id}'
            )
        if self.fm is None or self.fm.central_wl is None:
            cc.reset_energy(self.chunk)
        else:
            self.build_chunk_energy()
        self._logger.debug('End update_energy')

    def update_position(self):
        # DB - 01-18-19 - NIFS has no WCS info in
        # extension; use primary header
        #
        # DB - 04-03-19 - NIFS spatial WCS info in the
        # header has way too large a FOV so hard-code
        # this to the instrument's tiny 3" x 3" FOV.

        self._logger.debug('Begin update_position')
        n_axis1 = self._headers[-1]['NAXIS1']
        n_axis2 = self._headers[-1]['NAXIS2']

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

        header = self._headers[0]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(header)
        header['CRVAL2'] = self.get_dec(header)
        header['CDELT1'] = RADIUS_LOOKUP[self._name]
        header['CDELT2'] = RADIUS_LOOKUP[self._name]
        header['CROTA1'] = 0.0
        # DB 05-03-19 - persist NAXIS values for NIFS
        header['NAXIS1'] = n_axis1
        header['NAXIS2'] = n_axis2
        # DB 18-05-21
        # EQUINOX and RADESYS values should be fine for all NIFS files since
        # that’s what is used for all other Gemini data as far as I’m aware.
        header['EQUINOX'] = 2000.0
        header['RADESYS'] = 'FK5'

        header['CRPIX1'] = InstrumentType.get_crpix1(header)
        header['CRPIX2'] = InstrumentType.get_crpix2(header)
        header['CD1_1'] = self.get_cd11(header)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(header)

        wcs_parser = fits2caom2.WcsParser(
            header, self.obs_id, self._extension
        )
        if self.chunk is None:
            self.chunk = Chunk()
        wcs_parser.augment_position(self.chunk)
        self.chunk.position_axis_1 = 1
        self.chunk.position_axis_2 = 2
        self._logger.debug('End update_position')


class Niri(InstrumentType):
    def __init__(self, name):
        super(Niri, self).__init__(name)

    def get_filter_name(self):
        # NIRI - prefer header keywords
        result = InstrumentType._search_through_keys(
            self._headers[0], ['FILTER'], FILTER_VALUES_TO_IGNORE,
        )
        filter_names = result
        self.filter_name = filter_names

    def reset_energy(self, observation_type):
        return (
            super(Niri, self).reset_energy(observation_type) or (
                self.filter_name is not None and (
                    self.filter_name == '' or 'INVALID' in self.filter_name
                )
            )
        )

    def update_energy(self):
        self._logger.debug('Begin update_energy')
        reset_energy = False

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

        # DB 22-08-19
        # Guesstimate from changes in values for L filter/grism centered and
        # blue combinations.
        # M: 1100 for f6-2pixBl, 860 for f6-4pixBl and 490 for f6-6pixBl.
        # Skip energy for INVALID f-ratio.

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
                'f32-9pix': 450.0,  # f32-10pix
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
                'f32-9pix': 500.0,  # f32-10pix
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
                'f6-6pix': 460.0,
                'f6-2pixBl': 1100.0,
                'f6-4pixBl': 860.0,
                'f6-6pixBl': 490.0,
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
                'f32-9pix': 570.0,  # f32-10pix
            },
        }

        # NIRI.Hgrismf32-G5228

        # No energy information is determined for darks.  The
        # latter are sometimes only identified by a 'blank' filter.  e.g.
        # NIRI 'flats' are sometimes obtained with the filter wheel blocked
        # off.

        # The focal_plane_mask has values like f6-2pixBl (2-pixel wide blue
        # slit used with f/6 camera) or f32-6pix (6-pixel wide slit [centered]
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

        if 'Jcon(112)_G0235' in self.filter_name:
            # DB - 01-04-19 The G0235 filter is listed as ‘damaged’ on the
            # Gemini NIRI filters web site:
            # https://www.gemini.edu/sciops/instruments/niri/imaging/filters.
            # Not enough info is given there for SVO to add this filter to
            # their system.  Hardcode a central wavelength of 1.1232 microns
            # and a FWHM of 0.0092 microns
            filter_md = svofps.FilterMetadata('NIRI')
            filter_md.central_wl = 1.1232
            filter_md.bandpass = 0.0092
        elif 'Msort' in self.filter_name or 'Mgrism' in self.filter_name:
            # DB - 15-04-19
            # The NIRI blocking filter page,
            # https://www.gemini.edu/sciops/instruments/niri/spectroscopy/blocking-
            # filters, gives no transmission data for this filter for the SVO
            # to use.  Will need to hard-code it to central wavelength of
            # (4.4+6)/2 and width of 6-4.4 microns.

            filter_md = svofps.FilterMetadata('NIRI')
            filter_md.set_central_wl(6.0, 4.4)
            filter_md.set_bandpass(6.0, 4.4)
        else:
            filter_md = svofps.get_filter_metadata(
                self._name,
                self.filter_name,
                json_lookup.get('telescope'),
            )
            if filter_md is None:
                raise mc.CadcException(
                    f'{self._name}: mystery filter {self.filter_name} for '
                    f'{self.obs_id}'
                )

        self.filter_name = json_lookup.get('filter_name')
        if self.data_product_type == DataProductType.IMAGE:
            self._logger.debug(f'SpectralWCS imaging for {self.obs_id}.')
            self.fm = filter_md
            self.fm.adjust_resolving_power()
        elif self.data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(f'SpectralWCS spectroscopy for {self.obs_id}.')
            self.fm = svofps.FilterMetadata('NIRI')
            self.fm.central_wl = filter_md.central_wl
            self.fm.bandpass = filter_md.bandpass
            # add the 'split' call because NIRI: Mystery disperser value
            # Mgrism_G5206 for GN-2007B-Q-75-61-003
            disperser = json_lookup.get('disperser').split('_')[0]
            if disperser in [
                'Jgrism',
                'Jgrismf32',
                'Hgrism',
                'Hgrismf32',
                'Kgrism',
                'Kgrismf32',
                'Lgrism',
                'Mgrism',
            ]:
                bandpass_name = disperser[0]
                f_ratio = json_lookup.get('focal_plane_mask')
                self._logger.debug(
                    f'Bandpass name is {bandpass_name} f_ratio is '
                    f'{f_ratio} for {self.obs_id}'
                )
                if (
                        bandpass_name in NIRI_RESOLVING_POWER
                        and f_ratio in NIRI_RESOLVING_POWER[bandpass_name]
                ):
                    self.fm.resolving_power = NIRI_RESOLVING_POWER[bandpass_name][
                        f_ratio
                    ]
                elif 'pinhole' in f_ratio:
                    self._logger.info(
                        f'Pinhole. Setting energy to None for {self.obs_id}'
                    )
                    reset_energy = True
                elif f_ratio == 'INVALID':
                    self._logger.info(
                        f'INVALID f_ratio. Setting energy to None for '
                        f'{self.obs_id}'
                    )
                    reset_energy = True
                else:
                    raise mc.CadcException(
                        f'NIRI: Mystery bandpass name {bandpass_name} or '
                        f'f_ratio {f_ratio} for {self.obs_id}.'
                    )
            else:
                raise mc.CadcException(
                    f'NIRI: Mystery disperser value {disperser} for '
                    f'{self.obs_id}'
                )
        else:
            raise mc.CadcException(
                f'NIRI: Do not understand mode {self.data_product_type} for '
                f'{self.obs_id}'
            )

        if reset_energy:
            cc.reset_energy(self.chunk)
        else:
            self.build_chunk_energy()
        logging.debug('End _update_chunk_energy_niri')

    def update_position(self):
        self._logger.info(f'Begin update_position for {self.obs_id}')
        # DB 07-06-21
        # The extension CD values that are very, very close to 0 cause the
        # problems with Spatial WCS:
        # ERROR: spherepoly_from_array: a line segment overlaps or polygon too
        #        large
        # Try to use the primary values if this error occurs - there's an extra
        # '5' in the exponent
        if len(self._headers) > 1:
            pdu = self._headers[0]
            hdu0 = self._headers[1]
            pdu_cd1_1 = pdu.get('CD1_1')
            hdu0_cd1_1 = hdu0.get('CD1_1')
            if not math.isclose(pdu_cd1_1, hdu0_cd1_1) and math.isclose(
                    pdu_cd1_1 * 1e-50, hdu0_cd1_1
            ):
                pdu['NAXIS1'] = hdu0.get('NAXIS1')
                pdu['NAXIS2'] = hdu0.get('NAXIS2')
                wcs_parser = fits2caom2.WcsParser(
                    pdu, self.obs_id, self._extension
                )
                if self.chunk is None:
                    self.chunk = Chunk()
                wcs_parser.augment_position(self.chunk)
                if self.chunk.position is not None:
                    self.chunk.position_axis_1 = 1
                    self.chunk.position_axis_2 = 2
                    self.chunk.position.coordsys = pdu.get('FRAME')
                    self.chunk.position.equinox = mc.to_float(
                        pdu.get('EQUINOX')
                    )
        self._logger.info('End update_position')


class Oscir(InstrumentType):
    def __init__(self, name):
        super(Oscir, self).__init__(name)

    def get_art_product_type(self, header):
        # TODO - don't know what else to do right now
        return ProductType.SCIENCE

    def get_calibration_level(self, uri):
        result = CalibrationLevel.RAW_STANDARD
        if mc.CaomName(uri.lower()).file_id.startswith('r'):
            # DB 23-02-21
            # The best thing to do with OSCIR 'r' files is to add them as a
            # second cal level 2 plane and use the same metadata as the
            # unprocessed plane.
            result = CalibrationLevel.CALIBRATED
        return result

    def get_cd11(self, header):
        return RADIUS_LOOKUP[self._name]

    def get_cd22(self, header):
        return self.get_cd11(header)

    def get_dec(self, header):
        ra_ignore, dec = InstrumentType.get_sky_coord(
            header, 'RA_TEL', 'DEC_TEL'
        )
        return dec

    def get_exposure(self, header):
        # DB - 20-02-19 - json ‘exposure_time’ is in minutes, so multiply
        # by 60.0.
        return json_lookup.get('exposure_time') * 60.0

    def get_ra(self, header):
        ra, dec_ignore = InstrumentType.get_sky_coord(
            header, 'RA_TEL', 'DEC_TEL'
        )
        return ra

    def get_time_function_val(self, header):
        # DB - 20-02-19 - OSCIR json ‘ut_datetime’ is not correct.  Must
        # concatenate DATE-OBS and UTC1 values and convert to MJD as usual
        # or use MJD directly (seems to be correct starting MJD)
        #
        # DB - 06-03-19 json ut_date_time Gemini is based on DATE keyword.
        # Better to use MJD header keyword value directly.
        return header.get('MJD')

    def make_axes_consistent(self):
        # DB - 01-02-21
        # OSCIR data should all have NAXIS=6 in the header. Axis 1/2 are
        # position axes.  Axis 3 is the number of chop positions (1 or usually
        # = 2), axis 5 is the number of nod positions (1 or usually = 2), axis
        # 4 gives the number of ‘savesets’ per nod position and axis 6 gives
        # the number of ‘nod sets’.
        #
        # For that particular example you gave me I think the NAXIS6 = 1 is
        # incorrect.  Think it should equal the value of NODSETS.  And TOTFRMS
        # = 360 should equal naxis1 x naxis2 x naxis3 x naxis4.
        #
        # Spatial cutouts would use axes 1 and 2, although the precise
        # positions of each chop/nod position are not captured by the CAOM2
        # data.
        super(Oscir, self).make_axes_consistent()
        if self.chunk.naxis is not None and self.chunk.naxis == 6:
            self.chunk.naxis = 2
            self.chunk.time_axis = None

    def update_energy(self):
        """oscir-specific chunk-level Energy WCS construction."""
        self._logger.debug(f'Begin update_energy {self._name}')

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
        oscir_lookup = {
            'S_7.9': [7.91, 0.755],
            'S_8.8': [8.81, 0.871],
            'S_9.8': [9.80, 0.952],
            'S_10.3': [10.27, 1.103],
            'S_11.7': [11.70, 1.110],
            'S_12.5': [12.49, 1.156],
            'N_wide': [10.75, 5.230],
            'IHW_(17-19)': [18.17, 1.651],
            'Q3': [20.8, 1.650],
        }

        temp = self.filter_name
        self.filter_name = temp.split()[0]
        if self.filter_name not in oscir_lookup:
            raise mc.CadcException(
                f'{self._name}: Mystery FILTER keyword {self.filter_name} '
                f'for {self.obs_id}'
            )
        # DB 23-02-21
        # OSCIR files with the 'r' prefix will cause issues because they are
        # classified as spectra so remove the check
        #   'if data_product_type == DataProductType.IMAGE:'
        # and assume they are all images.  The filter bandpass is the best we
        # can do for the energy WCS in any case.
        self._logger.debug(f'oscir: SpectralWCS imaging mode for {self.obs_id}.')
        self.fm = svofps.FilterMetadata()
        self.fm.central_wl = oscir_lookup[self.filter_name][0]
        self.fm.bandpass = oscir_lookup[self.filter_name][1]
        self.build_chunk_energy()
        self._logger.debug(f'End update_energy {self._name}')

    def update_position(self):
        header = self._headers[self._extension]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(header)
        header['CRVAL2'] = self.get_dec(header)
        header['CDELT1'] = RADIUS_LOOKUP[self._name]
        header['CDELT2'] = RADIUS_LOOKUP[self._name]
        header['CROTA1'] = 0.0
        header['CRPIX1'] = Gpi.get_crpix1(header)
        header['CRPIX2'] = Gpi.get_crpix2(header)
        header['CD1_1'] = self.get_cd11(header)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(header)

        wcs_parser = fits2caom2.WcsParser(
            header, self.obs_id, self._extension
        )
        if self.chunk is None:
            self.chunk = Chunk()
        wcs_parser.augment_position(self.chunk)
        self.chunk.position_axis_1 = 1
        self.chunk.position_axis_2 = 2

        self.chunk.position.coordsys = header.get('FRAMEPA')
        self.chunk.position.equinox = mc.to_float(header.get('EQUINOX'))
        logging.debug('End _update_chunk_position')


class Phoenix(InstrumentType):
    def __init__(self, name):
        super(Phoenix, self).__init__(name)

    def get_art_product_type(self, header):
        result = ProductType.SCIENCE
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
            result = ProductType.CALIBRATION
        return result

    def get_calibration_level(self, uri):
        result = CalibrationLevel.RAW_STANDARD
        if mc.CaomName(uri.lower()).file_id.startswith('p'):
            result = CalibrationLevel.CALIBRATED
        return result

    def get_filter_name(self):
        super(Phoenix, self).get_filter_name()
        if self._filter_name is None or len(self._filter_name.strip()) == 0:
            result = InstrumentType._search_through_keys(
                self._headers[self._extension],
                ['CVF_POS', 'FILT_POS'],
                FILTER_VALUES_TO_IGNORE,
            )
            self._filter_name = result

    def get_obs_type(self, header):
        # DB - 18-02-19 - make PHOENIX searches more useful by making
        # estimated guesses at observation type
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
        object_value = json_lookup.get('object').lower()
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

    def reset_energy(self, observation_type):
        # DB 07-06-21
        # set energy to None
        return (
            super(Phoenix, self).reset_energy(observation_type) or
            (self.filter_name is not None and 'invalid' in self.filter_name)
        )

    def reset_position(self, headers, observation_type):
        result = super(Phoenix, self).reset_position(
            headers, observation_type
        )

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

        ra = json_lookup.get('ra')
        dec = json_lookup.get('dec')
        if ra is None and dec is None:
            result = True
        elif (
                ra is not None
                and math.isclose(ra, 0.0)
                and dec is not None
                and math.isclose(dec, 0.0)
        ):
            result = True
        return result

    def make_axes_consistent(self):
        super(Phoenix, self).make_axes_consistent()
        if self.chunk.naxis == 4:
            if self.chunk.position is None:
                self.chunk.naxis = None
                if self.chunk.time_axis is not None:
                    self.chunk.time_axis = None
            else:
                self.chunk.naxis = 2
                self.chunk.time_axis = None

        # DB - 05-06-20
        # That’s a composite observation (but with no way of determining the 4
        # members from the header) and it is an extracted spectrum and
        # (despite some header info suggesting otherwise) has no wavelength
        # scale. LINEAR must be a reference to the fact that the spacing of
        # each pixel is constant.  Since the scale is simply in pixels….
        # CTYPE1 will refer to the energy axis.   Would likely make more sense
        # for a value of PIXEL.
        ctype = self._headers[0].get('CTYPE1')
        if ctype is not None and ctype in [
            'LINEAR',
            'PIXEL',
        ]:
            self.chunk.naxis = None
            self.chunk.position_axis_1 = None
            self.chunk.position_axis_2 = None
            self.chunk.time_axis = None
            self.chunk.energy_axis = None

    def update_energy(self):
        self._logger.debug('Begin update_energy')
        # DB - 12-02-19 - Phoenix should be the same as TReCS but uses
        # NAXIS2 for the length of the dispersion axis.

        # select filter_id, wavelength_central, wavelength_lower,
        #     wavelength_upper
        # from gsa..gsa_filters where instrument = 'PHOENIX'
        # 0 - central
        # 1 - lower
        # 2 - upper
        #
        # units are microns
        #
        # dict with the filter wheels stripped from the names as returned by
        # query
        #
        # DB 30-04-19
        # A new filter not previously available.  Add info from here,
        # https://www.noao.edu/kpno/phoenix/filters.html

        PHOENIX = {
            '2030': [4.929000, 4.808000, 5.050000],
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
            '9440': [1.058500, 1.053000, 1.064000],
            'open': [3.0, 1.0, 5.0],
        }

        # DB - 12-02-19 - Note that the parenthetical numbers
        # after the Phoenix filter names (in the header) indicates the
        # filter wheel slot the filter is in and may occasionally change
        # so should be disregarded.
        if len(self.filter_name) > 0:
            self.filter_name = self.filter_name.split()[0]
            # found some files with '_' in the name
            if '_' in self.filter_name:
                self.filter_name = self.filter_name.split('_')[0]

        self._logger.debug(
            f'filter_name is {self.filter_name} for {self.obs_id}'
        )
        self.fm = svofps.FilterMetadata('Phoenix')
        if self.data_product_type in [
            DataProductType.SPECTRUM, DataProductType.IMAGE
        ]:
            self._logger.debug(
                f'DataProductType {self.data_product_type} for {self.obs_id}.'
            )
            if self.filter_name in PHOENIX:
                self.fm.set_bandpass(
                    PHOENIX[self.filter_name][2], PHOENIX[self.filter_name][1]
                )
                self.fm.central_wl = PHOENIX[self.filter_name][0]
            elif len(self.filter_name) == 0:
                # DB 11-02-21
                # With open filter in Phoenix the band pass coverage should
                # be from 1 to 5 microns, so central wavelength of 3 microns
                # and bandpass of 4 microns. Lines 2777 to 2778 should be
                # changed as well as adding an ‘open_(1)’ filter
                self.fm.set_bandpass(5.0, 1.0)
                self.fm.set_central_wl(5.0, 1.0)
            else:
                raise mc.CadcException(
                    f'mystery filter name {self.filter_name} for {self.obs_id}'
                )
        else:
            raise mc.CadcException(
                f'Unsupported DataProductType {self.data_product_type} for '
                f'{self.obs_id}'
            )
        self.build_chunk_energy()
        self._logger.debug('End update_energy')

    def update_position(self):
        self._logger.debug('Begin update_position')
        header = self._headers[self._extension]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(header)
        header['CRVAL2'] = self.get_dec(header)
        header['CDELT1'] = RADIUS_LOOKUP[self._name]
        header['CDELT2'] = RADIUS_LOOKUP[self._name]
        header['CROTA1'] = 0.0
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
        header['CRPIX1'] = InstrumentType.get_crpix1(header)
        header['CRPIX2'] = InstrumentType.get_crpix2(header)
        header['CD1_1'] = self.get_cd11(header)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(header)
        temp = header.get('EQUINOX')
        if temp is None or math.isclose(temp, 0.0):
            header['EQUINOX'] = header.get('EPOCH')

        wcs_parser = fits2caom2.WcsParser(
            header, self.obs_id, self._extension
        )
        if self.chunk is None:
            self.chunk = Chunk()
        wcs_parser.augment_position(self.chunk)
        self.chunk.position_axis_1 = 1
        self.chunk.position_axis_2 = 2
        self._logger.debug('End update_position')


class Texes(InstrumentType):
    def __init__(self, name):
        super(Texes, self).__init__(name)

    def get_calibration_level(self, uri):
        result = CalibrationLevel.RAW_STANDARD
        if '_red' in uri.lower() or '_sum' in uri.lower():
            result = CalibrationLevel.CALIBRATED
        return result

    def get_cd11(self, header, keyword='CDELT1'):
        return RADIUS_LOOKUP[self._name]

    def get_cd22(self, header):
        return self.get_cd11(header)

    def get_data_product_type(self, header):
        return DataProductType.SPECTRUM

    def get_dec(self, header):
        # bHROS, TEXES: ra/dec not in json
        return header.get('DEC')

    def get_exposure(self, header):
        # DB - 07-03-19 -  Use header ‘OBSTIME’ value as exposure time in
        # seconds.
        return mc.to_float(header.get('OBSTIME'))

    def get_filter_name(self):
        self._filter_name = None

    def get_ra(self, header):
        # bHROS, TEXES: ra/dec not in json
        return header.get('RA')

    def get_target_type(self):
        return TargetType.OBJECT

    def update_energy(self):
        # DB - 07-03-19
        # TEXES Spectroscopy
        #
        # Some special code will be needed for datalabels/planes.  There are
        # no datalabels in the FITS header.  json metadata (limited) must be
        # obtained with URL like
        # https://archive.gemini.edu/jsonsummary/canonical/
        #     filepre=TX20170321_flt.2507.fits.
        # Use TX20170321_flt.2507 as datalabel.  But NOTE:  *raw.2507.fits and
        # *red.2507.fits are two planes of the same observation. I’d suggest
        # we use ‘*raw*’ as the datalabel and ‘*red*’ or ‘*raw*’ as the
        # appropriate product ID’s for the science observations.  The ‘flt’
        # observations do not have a ‘red’ plane.  The json document contains
        # ‘filename’ if that’s helpful at all.  The ‘red’ files do not exist
        # for all ‘raw’ files.
        #
        #
        # Header OBSTYPE appears to be correct; not json obs_type.
        #
        # No previews are generated by Gemini
        #
        # Energy WCS:
        #
        # Central wavelength is given by 10,000/header(WAVENO0).  I have to do
        # some more investigation to see if we can determine wavelength
        # coverage (i.e. see if I can identify the echelle/echelon info from
        # the header - I don’t think so).  For now use 0.25 microns as the
        # fixed FWHM bandpass.
        self._logger.debug('Begin update_energy')
        if self.data_product_type == DataProductType.SPECTRUM:
            self._logger.debug(
                f'SpectralWCS spectral mode for {self.obs_id}.'
            )
            self.fm = svofps.FilterMetadata('TEXES')
            self.fm.central_wl = 10000 / self._headers[0].get('WAVENO0')
            self.fm.bandpass = 0.25
        else:
            # data_type/observing mode is always spectroscopy
            raise mc.CadcException(
                f'{self._name}: mystery data product type '
                f'{self.data_product_type} for {self.obs_id}'
            )
        self.build_chunk_energy()
        self._logger.debug('End update_energy')

    def update_position(self):
        header = self._headers[0]
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = self.get_ra(header)
        header['CRVAL2'] = self.get_dec(header)
        header['CDELT1'] = RADIUS_LOOKUP[self._name]
        header['CDELT2'] = RADIUS_LOOKUP[self._name]
        header['CROTA1'] = 0.0
        header['NAXIS1'] = 1
        header['NAXIS2'] = 1
        header['CRPIX1'] = InstrumentType.get_crpix1(header)
        header['CRPIX2'] = InstrumentType.get_crpix2(header)
        header['CD1_1'] = self.get_cd11(header)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CD2_2'] = self.get_cd22(header)
        wcs_parser = fits2caom2.WcsParser(
            header, self.obs_id, self._extension)
        if self.chunk is None:
            self.chunk = Chunk()
        wcs_parser.augment_position(self.chunk)
        self.chunk.position_axis_1 = 1
        self.chunk.position_axis_2 = 2
        self._logger.debug('End _update_chunk_position')


class Trecs(InstrumentType):
    def __init__(self, name):
        super(Trecs, self).__init__(name)

    def get_obs_intent(self, header):
        self._get_obs_class(header)
        result = ObservationIntentType.CALIBRATION
        if self._obs_class is None:
            self._get_obs_type(header)
            if self._obs_type in CAL_VALUES:
                result = ObservationIntentType.CALIBRATION
            else:
                result = ObservationIntentType.SCIENCE
                data_label = header.get('DATALAB')
                if data_label is not None and '-CAL' in data_label:
                    result = ObservationIntentType.CALIBRATION
        elif 'science' in self._obs_class:
            result = ObservationIntentType.SCIENCE
        return result

    def make_axes_consistent(self):
        super(Trecs, self).make_axes_consistent()
        if self.chunk.naxis == 4:
            if self.chunk.position is None:
                self.chunk.naxis = None
                if self.chunk.time_axis is not None:
                    self.chunk.time_axis = None
            else:
                self.chunk.naxis = 2
                self.chunk.time_axis = None

    def reset_energy(self, observation_type):
        # DB 23-04-19
        # GMOS GS-2005A-Q-26-12-001.  Lots of missing metadata, including
        # release date so no energy.  Ditto for TReCS GS-CAL20041206-6-007.
        return (
            super(Trecs, self).reset_energy(observation_type) or
            (self.filter_name is not None and self.filter_name == '')
        )

    def update_energy(self):
        self._logger.debug('Begin _update_chunk_energy_trecs')

        # Look at the json ‘disperser’ value.  If = LowRes-20" then
        # resolving power = 80.  If LowRes-10" then resolving power = 100.
        # There was a high-res mode but perhaps never used.  Again,
        # naxis = NAXIS1 header value and other wcs info as above for NIRI.
        # But might take some string manipulation to match filter names with
        # SVO filters.

        if self.filter_name == 'Datum+Datum':
            # DB 22-08-19
            # I think it might be better to not provide wavelength info rather
            # than a very wide range.  That’s how some other TReCS data with
            # null filter names are treated.
            self._logger.info(
                f'Filter is {self.filter_name}, no Spectral WCS for '
                f'{self.obs_id}'
            )
        else:
            if self.filter_name == 'Qone-17.8um':
                # DB 22-08-19
                # Qone-17.8um filter.  I can find no info about the bandpass
                # for that filter on the web but this info was in the old
                # gsa_filters table. Hardcode lower/upper bandpasses of 17.3
                # and 18.17 microns.
                w_min = 17.3
                w_max = 18.17
                self.fm = svofps.FilterMetadata()
                self.fm.set_bandpass(w_max, w_min)
                self.fm.set_central_wl(w_max, w_min)
                self.fm.set_resolving_power(w_max, w_min)
            else:
                self.fm = svofps.get_filter_metadata(
                    self._name,
                    self.filter_name,
                    json_lookup.get('telescope'),
                )
            if self.fm is None:
                raise mc.CadcException(
                    f'{Inst.TRECS}: Mystery filter '
                    f'{self.filter_name}'
                )
            if self.data_product_type == DataProductType.IMAGE:
                self._logger.debug(f'imaging mode for {self.obs_id}.')
            elif self.data_product_type == DataProductType.SPECTRUM:
                self._logger.debug(f'LS|Spectroscopy mode for {self.obs_id}.')
                disperser = json_lookup.get('disperser')
                if disperser is not None:
                    if disperser == 'LowRes-20':
                        self.fm.resolving_power = 80.0
                    elif disperser == 'LowRes-10':
                        self.fm.resolving_power = 100.0
            else:
                raise mc.CadcException(
                    f'TReCS: Do not understand mode {self.data_product_type} '
                    f'for {self.obs_id}'
                )
            self.build_chunk_energy()
        self._logger.debug('End _update_chunk_energy_trecs')

    def update_position(self):
        # DB 22-08-19
        # For the file with CTYPE1 = 0 in HDU 1, all
        # of the other WCS info appears to be fine
        # (i.e. identical with the primary header).
        # Use the HDU 0 values.
        if len(self._headers) > 1:
            ctype1 = self._headers[self._extension].get('CTYPE1')
            if isinstance(ctype1, str):
                # value repair for a small subset of TReCS files  :(
                # test is rS20060306S0090, GS-2005B-Q-10-63-003
                if ctype1 == '0':
                    self._headers[self._extension]['CTYPE1'] = 'RA---TAN'
                wcs_parser = fits2caom2.WcsParser(
                    self._headers[self._extension],
                    self.obs_id,
                    self._extension,
                )
            else:
                wcs_parser = fits2caom2.WcsParser(
                    self._headers[0], self.obs_id, self._extension
                )
        else:
            wcs_parser = fits2caom2.WcsParser(
                self._headers[0], self.obs_id, self._extension
            )

        if self.chunk is None:
            self.chunk = Chunk()

        wcs_parser.augment_position(self.chunk)
        if self.chunk.position is not None:
            self.chunk.position_axis_1 = 1
            self.chunk.position_axis_2 = 2
            self.chunk.position.coordsys = self._headers[0].get('FRAME')
            self.chunk.position.equinox = mc.to_float(
                self._headers[0].get('EQUINOX')
            )


def instrument_factory(name):
    lookup = {
        Inst.BHROS: Bhros,
        Inst.CIRPASS: Cirpass,
        Inst.F2: F2,
        Inst.FLAMINGOS: Flamingos,
        Inst.ALOPEKE: Fox,
        Inst.ZORRO: Fox,
        # Inst.GMOS: Gmos,
        # Inst.GMOSS: Gmos,
        # Inst.GMOSN: Gmos,
        Inst.GNIRS: Gnirs,
        Inst.GPI: Gpi,
        Inst.GRACES: Graces,
        Inst.GSAOI: Gsaoi,
        Inst.HOKUPAA: Hokupaa,
        Inst.HRWFS: Hrwfs,
        Inst.IGRINS: Igrins,
        Inst.MICHELLE: Michelle,
        Inst.NICI: Nici,
        Inst.NIFS: Nifs,
        Inst.NIRI: Niri,
        Inst.OSCIR: Oscir,
        Inst.PHOENIX: Phoenix,
        Inst.TEXES: Texes,
        Inst.TRECS: Trecs,
    }
    if name in lookup:
        return lookup.get(name)(name)
    else:
        raise mc.CadcException(f'Mystery name {name}.')
