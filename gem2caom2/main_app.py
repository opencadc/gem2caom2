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

import logging
import math
import re
import traceback

from caom2 import Observation, CalibrationLevel, Chunk, ProductType
from caom2 import TypedList, DerivedObservation, DataProductType
from caom2 import ObservationIntentType, TargetType
from caom2utils import WcsParser, update_artifact_meta
from caom2pipe import manage_composable as mc
from caom2pipe import caom_composable as cc
from caom2pipe import astro_composable as ac

import gem2caom2.external_metadata as em
import gem2caom2.obs_file_relationship as ofr
from gem2caom2.gem_name import GemName
from gem2caom2 import instruments, program_metadata, svofps
from gem2caom2.gemini_reader import GeminiMetadataLookup
from gem2caom2.util import Inst, COLLECTION, SCHEME


__all__ = ['GeminiMapping', 'APPLICATION', 'mapping_factory']

APPLICATION = 'gem2caom2'

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


class GeminiMapping(cc.TelescopeMapping):

    def __init__(self, storage_name, headers, lookup, instrument):
        super().__init__(storage_name, headers)
        self._metadata_reader = lookup.reader
        self._instrument = instrument
        self._lookup = lookup

    def get_time_delta(self, ext):
        exptime = self.get_exposure(ext)
        if exptime is None:
            return None
        return mc.to_float(exptime) / (24.0 * 3600.0)

    def get_calibration_level(self, ext):
        result = CalibrationLevel.RAW_STANDARD
        reduction = self._lookup.reduction(self._storage_name.file_uri)
        if (
            reduction is not None and (
                'PROCESSED' in reduction or 'PREPARED' in reduction
            )
        ):
            result = CalibrationLevel.CALIBRATED
        return result

    def get_art_product_type(self, ext):
        """
        If obsclass is unknown then CAOM2 ProductType is set to CALIBRATION.
        This should only effect early data when OBSCLASS was not in the JSON
        summary metadata or FITS headers.
        """
        data_label = self._lookup.data_label(self._storage_name.file_uri)
        obs_type = self._lookup.observation_type(self._storage_name.file_uri)
        obs_class = self._lookup.observation_class(self._storage_name.file_uri)

        self._logger.debug(
            f'obs type is {obs_type} obs class is {obs_class} for {data_label}'
        )
        if obs_type is not None and obs_type == 'MASK':
            result = ProductType.AUXILIARY
        elif obs_class is None:
            # 'unknown' is the value of observation_class for CIRPASS
            if (
                obs_type is not None and
                obs_type in ['OBJECT', 'unknown']
            ):
                if data_label is not None and 'CAL' in data_label:
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

    def get_data_product_type(self, ext):
        mode = self._lookup.mode(self._storage_name.file_uri)
        obs_type = self._lookup.observation_type(self._storage_name.file_uri)

        if mode is None:
            raise mc.CadcException(
                f'No mode information found for {self._storage_name.file_name}'
            )
        elif (mode == 'imaging') or (
            obs_type is not None and obs_type == 'MASK'
        ):
            result = DataProductType.IMAGE
        else:
            result = DataProductType.SPECTRUM
        return result

    def get_data_release(self, ext):
        """
        Determine the plane-level data release date.
        :return: The Plane release date, or None if not found.
        """
        # every instrument has a 'release' keyword in the JSON summary
        # not every instrument (Michelle) has a RELEASE keyword in
        # the appropriate headers
        result = self._lookup.release(self._storage_name.file_uri)
        if result is not None and result.startswith('0001'):
            # because obs id GN-2008A-Q-39-69-015
            result = result.replace('0001', '2001')
        return result

    def get_dec(self, ext):
        return em.current_instrument.get_dec(self._headers[ext])

    def get_exposure(self, ext):
        """
        Calculate the exposure time. EXPTIME in the header is not always the
        total exposure time for some of the IR instruments.  For these EXPTIME
        is usually the individual exposure time for each chop/nod and a bunch
        of chops/nods are executed for one observation.  Gemini correctly
         allows for this in the json ‘exposure_time’ value.
        """
        return self._lookup.exposure_time(self._storage_name.file_uri)

    # def get_instrument_name(self, ext):
    #     return em.repair_instrument(
    #         self._lookup.instrument(self._storage_name.file_uri)
    #     ).value

    def get_meta_release(self, ext):
        """
        Determine the metadata release date (Observation and Plane-level).
        :return: The Observation/Plane release date, or None if not found.
        """
        # TODO - this is a bit different ROTFL
        if self._headers is None:
            # GenericParser, so no headers retrieved from archive.gemini.edu,
            # probably a 403 being returned by the site, assume proprietary
            meta_release = self._lookup.release(self._storage_name.file_uri)
        else:
            # DB 21-08-19
            # If PROP_MD is T, use JSON ‘release’ value for metadata release
            # date. If no PROP_MD present or value is F use the JSON
            # ut_datetime value.
            prop_md = self._headers[ext].get('PROP_MD')
            if prop_md is None or prop_md is False or prop_md == 'F':
                meta_release = self._lookup.ut_datetime(
                    self._storage_name.file_uri
                )
            else:
                meta_release = self._lookup.release(self._storage_name.file_uri)
        return meta_release

    def get_obs_intent(self, ext):
        result = ObservationIntentType.CALIBRATION
        self._get_data_label()
        self._get_obs_class(ext)
        self._logger.debug(
            f'observation_class is {self._obs_class} for {self._data_label}'
        )
        if self._obs_class is None:
            obs_type = self._lookup.observation_type(
                self._storage_name.file_uri
            )

            self._logger.debug(
                f'observation_type is {self._obs_type} for {self._data_label}'
            )
            if obs_type is None:
                self._logger.debug(f'data_label is {self._data_label}')
                if (
                    self._data_label is None or (
                        self._data_label is not None and
                        '-CAL' not in self._data_label
                    )
                    and self._headers[ext] is not None
                ):
                    object_value = self._headers[ext].get('OBJECT')
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
                if obs_type in CAL_VALUES:
                    result = ObservationIntentType.CALIBRATION
                else:
                    result = ObservationIntentType.SCIENCE
        elif 'science' in self._obs_class:
            result = ObservationIntentType.SCIENCE
        return result

    def get_obs_type(self, ext):
        result = self._lookup.observation_type(self._storage_name.file_uri)
        obs_class = self._lookup.observation_class(self._storage_name.file_uri)
        if (
            obs_class is not None and
            (obs_class == 'acq' or obs_class == 'acqCal')
        ):
            result = 'ACQUISITION'
        return result

    def get_proposal_id(self, ext):
        return self._lookup.program_id(self._storage_name.file_uri)

    def get_provenance_keywords(self, ext):
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
        return self._lookup.mode(self._storage_name.file_uri)

    def get_provenance_last_executed(self, ext):
        def breakout(comments):
            result = None
            temp = comments.split('\n')
            if len(temp) > 6 and 'HST' in temp[6]:
                # go from HST to UTC
                result = mc.make_time(temp[6])
            return result

        return self._get_provenance_breakout(ext, breakout)

    def get_provenance_producer(self, ext):
        cal_level = self.get_calibration_level(self._storage_name.file_uri)
        if cal_level in [
            CalibrationLevel.CALIBRATED,
            CalibrationLevel.PRODUCT,
            CalibrationLevel.ANALYSIS_PRODUCT,
        ]:
            comments = str(self._headers[ext].get('COMMENT'))
            result = comments.split('Processed by the')[1].split('|')[0]
        else:
            result = self._headers[ext].get('IMAGESWV')
        return result

    def get_provenance_reference(self, ext):
        def breakout(comments):
            return 'https://www.gemini.edu/instrumentation/graces/data-reduction'

        return self._get_provenance_breakout(ext, breakout)

    def get_provenance_version(self, ext):
        def breakout(comments):
            temp = comments.split('opera-')[1].split('build date')[0]
            return f'opera-{temp}'

        return self._get_provenance_breakout(ext, breakout)

    def _get_provenance_breakout(self, ext, fn):
        result = None
        cal_level = self.get_calibration_level()
        if cal_level in [
            CalibrationLevel.CALIBRATED,
            CalibrationLevel.PRODUCT,
            CalibrationLevel.ANALYSIS_PRODUCT,
        ]:
            comments = str(self._headers[ext].get('COMMENT'))
            result = fn(comments)
        return result

    def get_ra(self, header):
        """
        Get the right ascension. Rely on the JSON metadata, because it's all in
        the same units (degrees).

        :param header:  The FITS header for the current extension.
        :return: ra, or None if not found.
        """
        return self._lookup.ra(self._storage_name.file_uri)

    def get_target_moving(self, ext):
        """
        Calculate whether the Target moving.
        Non-sidereal tracking -> setting moving target to "True"

        DB, 01-08-19
        Many calibration observations are acquired with the telescope parked and
        hence not tracking at sidereal rate.

        :param header:  The FITS header for the current extension.
        :return: The Target TargetType, or None if not found.
        """
        types = self._lookup.types(self._storage_name.file_uri)
        if 'NON_SIDEREAL' in types:
            return True
        else:
            return None

    def get_target_type(self):
        """
        Calculate the Target TargetType
        """
        result = TargetType.FIELD
        spectroscopy = self._lookup.spectroscopy(self._storage_name.file_uri)
        if spectroscopy:
            result = TargetType.OBJECT
        return result

    def get_time_function_val(self, header):
        """
        Calculate the Chunk Time WCS function value, in 'mjd'.

        :param header:  The FITS header for the current extension (not used).
        :return: The Time WCS value from JSON Summary Metadata.
        """
        time_string = self._lookup.ut_datetime(self._storage_name.file_uri)
        return ac.get_datetime(time_string)

    def _get_data_label(self):
        return self._lookup.data_label(self._storage_name.file_uri)

    def accumulate_blueprint(self, bp, application=None):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp, APPLICATION)
        self._logger.debug(
            f'Begin accumulate_fits_bp for {self._storage_name.file_id}.'
        )
        em.defining_metadata_finder.get_obs_metadata(
            self._storage_name.file_id
        )

        bp.set('Observation.type', 'get_obs_type()')
        bp.set('Observation.intent', 'get_obs_intent()')
        bp.set('Observation.metaRelease', 'get_meta_release()')
        bp.set('Observation.target.type', 'get_target_type()')
        bp.set('Observation.target.moving', 'get_target_moving()')
        bp.set('Observation.proposal.id', 'get_proposal_id()')

        bp.clear('Observation.algorithm.name')
        bp.set('Observation.instrument.name', self._instrument.value)
        # instrument = em.get_instrument(self._storage_name.file_uri)
        logging.error(self._metadata_reader.json_metadata)
        telescope = self._lookup.telescope(self._storage_name.file_uri)
        logging.error('after')
        if telescope is not None:
            if telescope is not None and 'North' in telescope:
                x, y, z = ac.get_location(19.823806, -155.46906, 4213.0)
            else:
                x, y, z = ac.get_location(-30.240750, -70.736693, 2722.0)
            bp.set('Observation.telescope.geoLocationX', x)
            bp.set('Observation.telescope.geoLocationY', y)
            bp.set('Observation.telescope.geoLocationZ', z)

        bp.set('Plane.productID', self._storage_name.file_id)
        bp.set('Plane.dataProductType', 'get_data_product_type()')
        bp.set('Plane.calibrationLevel', 'get_calibration_level()')
        bp.set('Plane.metaRelease', 'get_meta_release()')
        bp.set('Plane.dataRelease', 'get_data_release()')

        bp.set('Plane.provenance.name', 'Gemini Observatory Data')
        bp.set('Plane.provenance.project', 'Gemini Archive')
        # Add IMAGESWV for GRACES
        bp.add_fits_attribute('Plane.provenance.producer', 'IMAGESWV')
        bp.set_default('Plane.provenance.producer', 'Gemini Observatory')
        data_label = self._lookup.data_label(self._storage_name.file_uri)
        bp.set(
            'Plane.provenance.reference',
            f'http://archive.gemini.edu/searchform/{data_label}',
        )

        # if instrument is Inst.GRACES:
        #     bp.set(
        #         'Plane.provenance.lastExecuted',
        #         'get_provenance_last_executed(parameters)',
        #     )
        #     bp.set(
        #         'Plane.provenance.producer',
        #         'get_provenance_producer(parameters)',
        #     )
        #     bp.set(
        #         'Plane.provenance.reference',
        #         'get_provenance_reference(parameters)',
        #     )
        #     bp.set(
        #         'Plane.provenance.version',
        #         'get_provenance_version(parameters)',
        #     )

        bp.set('Artifact.productType', 'get_art_product_type()')
        # bp.set('Artifact.contentChecksum', f'md5:{json_lookup.get("data_md5")}')
        # bp.set('Artifact.contentLength', json_lookup.get('data_size'))
        # bp.set('Artifact.contentType', 'application/fits')
        # always see the metadata, see the data only when it's public
        bp.set('Artifact.releaseType', 'data')
        bp.set('Artifact.uri', self._storage_name.file_uri)

        # # if instrument is Inst.CIRPASS:
        # #     bp.set_default('Observation.telescope.name', 'Gemini-South')
        # mode = json_lookup.get('mode')
        # if not (
        #     instrument
        #     in [
        #         Inst.GPI,
        #         Inst.PHOENIX,
        #         Inst.HOKUPAA,
        #         Inst.OSCIR,
        #         Inst.BHROS,
        #         Inst.TRECS,
        #     ]
        #     or (
        #         instrument is Inst.GRACES
        #         and mode is not None
        #         and mode != 'imaging'
        #     )
        # ):
        #     bp.configure_position_axes((1, 2))

        # if instrument is Inst.FLAMINGOS:
        #     # DB 27-05-19
        #     # Flamingos, you actually want to use the EQUINOX value, not the
        #     # EPOCH.   And I think EQUINOX header value is usually 2000.0, even
        #     # for the example GS-CAL20020620-15-0462 02jun20.0462 with
        #     # RA_TEL = “UNAVAILABLE”.  For Gemini the assumption is that the
        #     # RA/Dec values in the headers are always based on the position of
        #     # the equinox given at the time specified by the EQUINOX keyword
        #     # value.
        #     bp.clear('Chunk.position.equinox')
        #     bp.add_fits_attribute('Chunk.position.equinox', 'EQUINOX')

        bp.configure_time_axis(3)

        # The Chunk time metadata is calculated using keywords from the
        # primary header, and the only I could figure out to access keywords
        # in the primary is through a function. JB
        bp.set('Chunk.time.resolution', 'get_exposure()')
        bp.set('Chunk.time.exposure', 'get_exposure()')

        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.error.syser', '1e-07')
        bp.set('Chunk.time.axis.error.rnder', '1e-07')
        bp.set('Chunk.time.axis.function.naxis', '1')
        # if instrument in [Inst.ALOPEKE, Inst.ZORRO]:
        #     bp.clear('Chunk.time.axis.function.naxis')
        #     bp.add_fits_attribute('Chunk.time.axis.function.naxis', 'NAXIS3')
        #     bp.set_default('Chunk.time.axis.function.naxis', 1)

        bp.set('Chunk.time.axis.function.delta', 'get_time_delta()')
        bp.set('Chunk.time.axis.function.refCoord.pix', '0.5')
        bp.set(
            'Chunk.time.axis.function.refCoord.val',
            'get_time_function_val(header)',
        )

        self._logger.debug('Done accumulate_fits_bp.')

    def update(self, observation, file_info, caom_repo_client=None):
        """Called to fill multiple CAOM model elements and/or attributes, must
        have this signature for import_module loading and execution.

        :param observation A CAOM Observation model instance.
        :param **kwargs Everything else."""
        self._logger.debug('Begin update.')
        mc.check_param(observation, Observation)

        # processed files
        if cc.is_composite(self._headers) and not isinstance(
            observation, DerivedObservation
        ):
            observation = self._update_composite(observation)

        if self._instrument in [Inst.MICHELLE, Inst.GNIRS]:
            # DB 16-04-19
            # The more important issue with this and other files is that they
            # contain no image extensions.  The file is downloadable from
            # the Gemini archive but their only content is the primary
            # header.   There is no pixel data.  Test for the existence of a
            # FITS extension and skip processing of a michelle file if
            # there isn’t one

            # DB 18-04-19
            #
            # For the last GNIRS file (NAXIS=0)  skip the file if it doesn’t
            # have an extension.

            if len(self._headers) == 1:
                self._logger.warning(
                    f'{self._instrument}: no image data for '
                    f'{observation.observation_id}. Cannot build an '
                    f'observation.'
                )
                return None

        config = mc.Config()
        config.get_executors()
        try:
            for plane in observation.planes.values():
                if (
                    self._storage_name.product_id is not None
                    and self._storage_name.product_id != plane.product_id
                ):
                    continue

                for artifact in plane.artifacts.values():
                    self._should_artifact_be_renamed(artifact)
                    if GemName.is_preview(artifact.uri):
                        continue

                    update_artifact_meta(artifact, file_info)
                    caom_name = mc.CaomName(artifact.uri)
                    file_id = ofr.remove_extensions(
                        mc.CaomName(caom_name.uri).file_name
                    )
        #             self._metadata_reader.json_metadata.get(
        #     self._storage_name.file_uri
        # ).reset_index(file_id)
                    processed = ofr.is_processed(caom_name.file_name)
                    if self._instrument in [
                        Inst.MICHELLE,
                        Inst.TRECS,
                        Inst.GNIRS,
                    ]:
                        # Michelle is a retired visitor instrument.
                        # Spatial WCS info is in primary header. There
                        # are a variable number of FITS extensions
                        # defined by primary keyword NUMEXT; assume the
                        # same spatial WCS for each for now - it differs
                        # only slightly because of telescope 'chopping'
                        # and 'nodding' used in acquisition. DB - 01-18-19
                        #
                        # DB - 01-18-19 - GNIRS has no WCS info in extension;
                        # use primary header
                        self._update_position_from_zeroth_header(
                            artifact, observation.observation_id,
                        )

                    delete_these_parts = []
                    for part in artifact.parts:

                        if part == '2' and self._instrument is Inst.GPI:
                            # GPI data sets have two extensions. First is
                            # science image (with WCS), second is data quality
                            # for each pixel (no WCS).
                            self._logger.info(
                                f'GPI: Setting chunks to None for part {part} '
                                f'for {observation.observation_id}'
                            )
                            artifact.parts[part].chunks = TypedList(
                                Chunk,
                            )
                            continue
                        for c in artifact.parts[part].chunks:
                            # example is CIRPASS/2003jun30_3385.fits - older
                            # versions of the file had more headers
                            if int(part) >= len(self._headers):
                                delete_these_parts.append(part)
                                continue

                            # energy WCS
                            x = instruments.instrument_factory(self._instrument)
                            x.headers = self._headers
                            x.extension = int(part)
                            x.get_filter_name()
                            x.chunk = c
                            if x.reset_energy(observation.type):
                                cc.reset_energy(c)
                            else:
                                x.data_product_type = plane.data_product_type
                                x.obs_id = observation.observation_id
                                x.update_energy()

                            # position WCS
                            mode = self._lookup.mode(
                                self._storage_name.file_uri
                            )
                            x.mode = mode
                            if x.reset_position(self._headers, observation.type):
                                self._logger.debug(
                                    f'Setting Spatial WCS to None for '
                                    f'{observation.observation_id}'
                                )
                                cc.reset_position(c)
                            else:
                                x.update_position()

                            # time WCS
                            x.update_time()
                            x.make_axes_consistent()

                    if isinstance(observation, DerivedObservation):
                        values = cc.find_keywords_in_headers(
                            self._headers[1:], ['IMCMB']
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
                        or self._instrument is Inst.TEXES
                    ) and 'jpg' not in caom_name.file_name:
                        # not the preview artifact
                        if plane.provenance is not None:
                            if self._instrument is not Inst.GRACES:
                                plane.provenance.reference = (
                                    f'http://archive.gemini.edu/searchform/'
                                    f'filepre={caom_name.file_name}'
                                )

                    for part in delete_these_parts:
                        self._logger.warning(
                            f'Delete part {part} from artifact {artifact.uri}'
                        )
                        artifact.parts.pop(part)

                program = program_metadata.get_pi_metadata(
                    observation.proposal.id
                )
                if program is not None:
                    observation.proposal.pi_name = program['pi_name']
                    observation.proposal.title = program['title']

            if isinstance(observation, DerivedObservation):
                cc.update_observation_members(observation)

            em.value_repair.repair(observation)
        except Exception as e:
            self._logger.error(
                f'Error {e} for {observation.observation_id} instrument '
                f'{self._instrument}'
            )
            tb = traceback.format_exc()
            self._logger.debug(tb)
            raise mc.CadcException(e)
        self._logger.debug('Done update.')
        return observation

    def _should_artifact_be_renamed(self, artifact):
        if artifact.uri.startswith('gemini'):
            if artifact.uri.startswith('gemini:GEM/'):
                artifact.uri = artifact.uri.replace(
                    'gemini:GEM/', 'gemini:GEMINI/'
                )
        if artifact.uri.startswith('ad'):
            artifact.uri = artifact.uri.replace(
                'ad:GEM/', 'cadc:GEMINI/'
            )

    def _update_position_from_zeroth_header(self, artifact, obs_id):
        """Make the 0th header spatial WCS the WCS for all the
        chunks."""
        primary_header = self._headers[0]

        # naxis values are only available from extensions

        primary_header['NAXIS1'] = self._headers[-1].get('NAXIS1')
        primary_header['NAXIS2'] = self._headers[-1].get('NAXIS2')

        primary_chunk = Chunk()
        types = self._lookup.types(self._storage_name.file_uri)
        ra = self._headers[0].get('RA')
        ra_json = self.get_ra(0)
        dec_json = self._lookup.dec(self._storage_name.file_uri)
        if ('AZEL_TARGET' in types and ra is None) or (
            self._instrument is Inst.GNIRS and ra_json is None and
            dec_json is None
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

            self._logger.info(
                f'{self._instrument}: Spatial WCS is None for {obs_id}'
            )
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

    def _update_composite(self, obs):
        result = None
        if self._instrument is Inst.TRECS:
            if self._storage_name.product_id is not None and (
                self._storage_name.product_id.startswith('rS')
                or self._storage_name.product_id.startswith('rN')
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
            self._logger.info(f'{obs.observation_id} is a Derived Observation.')
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
    logging.debug(f'Being _repair_provenance_value of {value} for {obs_id}.')
    prov_file_id = value
    try:
        uri = mc.build_uri(COLLECTION, prov_file_id, SCHEME)
        temp = em.defining_metadata_finder.get(uri)
        prov_obs_id = temp.data_label
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


class Cirpass(GeminiMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def accumulate_blueprint(self, bp, application=None):
        super().accumulate_blueprint(bp)
        bp.set_default('Observation.telescope.name', 'Gemini-South')


class Flamingos(GeminiMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def accumulate_blueprint(self, bp, application=None):
        super().accumulate_blueprint(bp)
        # DB 27-05-19
        # Flamingos, you actually want to use the EQUINOX value, not the
        # EPOCH.   And I think EQUINOX header value is usually 2000.0, even
        # for the example GS-CAL20020620-15-0462 02jun20.0462 with
        # RA_TEL = “UNAVAILABLE”.  For Gemini the assumption is that the
        # RA/Dec values in the headers are always based on the position of
        # the equinox given at the time specified by the EQUINOX keyword
        # value.
        bp.clear('Chunk.position.equinox')
        bp.add_fits_attribute('Chunk.position.equinox', 'EQUINOX')


class Fox(GeminiMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def accumulate_blueprint(self, bp, application=None):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        bp.set(
            'Plane.provenance.reference',
            f'http://archive.gemini.edu/searchform/filepre='
            f'{self._storage_name.file_id}.fits',
        )
        bp.clear('Chunk.time.axis.function.naxis')
        bp.add_fits_attribute('Chunk.time.axis.function.naxis', 'NAXIS3')
        bp.set_default('Chunk.time.axis.function.naxis', 1)


class Gmos(GeminiMapping):

    def __init__(self, storage_name, headers, lookup, instrument):
        super().__init__(storage_name, headers, lookup, instrument)

    def accumulate_blueprint(self, bp, application=None):
        super().accumulate_blueprint(bp)
        bp.set(
            'Observation.instrument.keywords',
            'get_provenance_keywords(uri)',
        )

    def get_exposure(self, ext):
        if self.get_obs_type(ext) == 'MASK':
            result = 0.0
        else:
            result = super(Gmos, self).get_exposure(ext)
        return result

    def get_time_delta(self, ext):
        if self.get_obs_type(ext) == 'MASK':
            # DB - 05-03-19 - delta hardcoded to 0
            exptime = 0.0
        else:
            exptime = super(Gmos, self).get_time_delta(ext)
        return exptime

    def reset_position(self, headers, observation_type):
        # DB - 04-03-19
        # Another type of GMOS-N/S dataset to archive.
        # Mask images.   json observation_type = “MASK”.
        # These have no WCS info at all, although I guess
        # json ut_date_time could be used as the start date
        # with null exposure time. These would have only
        # instrument, obstype, datatype (spectrum) and
        # product type (AUXILIARY) set.
        return (
            super(Gmos, self).reset_position(headers, observation_type) or
            observation_type == 'MASK'
        )

    def update_energy(self):
        self._logger.debug('Begin _update_chunk_energy_gmos')

        GMOS_RESOLVING_POWER = {
            'B1200': 3744.0,
            'R831': 4396.0,
            'B600': 1688.0,
            'R600': 3744.0,
            'R400': 1918.0,
            'R150': 631.0,
        }
        # DB 02-12-19
        # Found basic info for the Lya395 filter in a published paper.
        # Central wavelength is .3955 microns and FWHM is 0.00327 microns.
        #
        # The Hartmann ‘filter’ is actually a mask that blocks off half of the
        # light path to the detector but is ‘open’ on the other half.  So same
        # ‘bandpass’ as ‘open’ or ‘empty’.    It’s used when focusing the camera.
        #
        # Blocking filters are only used when doing multi-object spectroscopy.  By
        # reducing the wavelength coverage they are able to add a second row of
        # slits to the mask since the spectra then don’t cover the entire width of
        # the detector in the dispersion direction.

        # 0 = min
        # 1 = max
        # units are microns
        lookup = {
            'GG455': [0.46000, 1.10000],
            'OG515': [0.52000, 1.10000],
            'RG610': [0.61500, 1.10000],
            'RG780': [0.7800, 1.10000],
            'Lya395': [0.393865, 0.397135],
            'open': [0.35000, 1.10000],
            'open2-8': [0.35000, 1.10000],
            'empty': [0.35000, 1.10000],
            'HartmannA': [0.35000, 1.1000],
            'HartmannB': [0.35000, 1.1000],
            'HartmannC': [0.35000, 1.1000],
            'HartmannD': [0.35000, 1.1000],
        }

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

        # DB 02-12-19
        #
        # One problem is 186 files with filter values of empty_0[134].  The
        # ‘empty*’ component should be treated like ‘open’ as far as bandpass is
        # concerned.   And ‘empty_*’ could be stripped out of the displayed filter
        # name. Likely less confusing to remove ‘empty*’ completely.

        reset_energy = False

        filter_md = None
        if (
            'open' not in self.filter_name
            and 'No_Value' not in self.filter_name
            and 'empty' not in self.filter_name
        ):
            filter_md = svofps.get_filter_metadata(
                self._name,
                self.filter_name,
                self._lookup.telescope(self._storage_name.file_uri),
            )

        if 'empty' in self.filter_name:
            # set to 'empty' string here, so can still use lookup logic
            self.filter_name = re.sub('empty_\\d*', 'empty', self.filter_name)

        w_max = 10.0
        w_min = 0.0
        for ii in self.filter_name.split('+'):
            if ii in lookup:
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

        filter_md = svofps.FilterMetadata()
        filter_md.set_bandpass(w_max, w_min)
        filter_md.set_central_wl(w_max, w_min)
        filter_md.set_resolving_power(w_max, w_min)

        if self.data_product_type == DataProductType.SPECTRUM:
            logging.debug(
                f'{self._name}: SpectralWCS spectroscopy for {self.obs_id}.'
            )
            if math.isclose(filter_md.central_wl, 0.0):
                logging.info(
                    f'{self._name}: no spectral wcs, central wavelength is '
                    f'{filter_md.central_wl} for {self.obs_id}'
                )
                return
            disperser = self._lookup.disperser(self._storage_name.file_uri)
            # 'unknown' in disperser test obs is GN-2004B-Q-30-15-002
            # 'OLDMIRROR' in disperser test obs is GN-CAL20020329-2-025
            # 'No Value' in disperser test obs is GN-2013B-SV-152-120-001
            if (
                disperser is None
                or 'unknown' in disperser
                or 'OLDMIRROR' in disperser
                or 'No Value' in disperser
            ):
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
                logging.warning(
                    f'{self._name}: disperser is {disperser}, no energy.'
                )
                reset_energy = True
            else:
                self.fm = svofps.FilterMetadata()
                self.fm.central_wl = filter_md.central_wl
                self.fm.bandpass = filter_md.bandpass
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
                    self.fm.resolving_power = GMOS_RESOLVING_POWER[disperser]
                else:
                    raise mc.CadcException(
                        f'{self._name}: mystery disperser {disperser} for '
                        f'{self.obs_id}'
                    )
        elif self.data_product_type == DataProductType.IMAGE:
            self._logger.debug(
                f'{self._name}: SpectralWCS imaging for {self.obs_id}.'
            )
            self.fm = filter_md
        else:
            raise mc.CadcException(
                f'{self._name}: mystery data product type '
                f'{self.data_product_type} for {self.obs_id}'
            )
        if reset_energy:
            cc.reset_energy(self.chunk)
        else:
            self.build_chunk_energy()
        self._logger.debug('End _update_chunk_energy_gmos')


class Graces(GeminiMapping):

    def __init__(self, storage_name, headers, metadata_reader):
        super().__init__(storage_name, headers, metadata_reader)

    def accumulate_blueprint(self, bp, application=None):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
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
        mode = self._lookup.mode(self._storage_name.file_uri)
        if mode is not None and mode != 'imaging':
            bp.configure_position_axes((1, 2))


class Igrins(GeminiMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def accumulate_blueprint(self, bp, application=None):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        telescope = self._metadata_reader.json_metadata.get(
            self._storage_name.file_uri
        ).get('telescope')
        if telescope is not None and 'North' in telescope:
            x, y, z = ac.get_location(19.823806, -155.46906, 4213.0)
        else:
            x, y, z = ac.get_location(-30.240750, -70.736693, 2722.0)
        bp.set('Observation.telescope.geoLocationX', x)
        bp.set('Observation.telescope.geoLocationY', y)
        bp.set('Observation.telescope.geoLocationZ', z)


class Texes(GeminiMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def accumulate_blueprint(self, bp, application=None):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        bp.clear('Plane.provenance.reference')


def mapping_factory(storage_name, headers, metadata_reader):
    metadata_lookup = GeminiMetadataLookup(metadata_reader)
    inst = metadata_lookup.instrument(storage_name.file_uri)
    # if inst is Inst.TEXES:
    #     return Texes(storage_name, headers)
    # elif inst is Inst.IGRINS:
    #     return Igrins(storage_name, headers)
    # elif inst is Inst.GRACES:
    #     return Graces(storage_name, headers)
    # elif inst in [Inst.GMOS, Inst.GMOSN, Inst.GMOSS]:
    #     return Gmos(storage_name, headers, lookup, inst)
    # elif inst in [Inst.ALOPEKE, Inst.ZORRO]:
    #     return Fox(storage_name, headers)
    # elif inst is Inst.FLAMINGOS:
    #     return Flamingos(storage_name, headers)
    # elif inst is Inst.CIRPASS:
    #     return Cirpass(storage_name, headers)
    # else:
    #     return GeminiMapping(storage_name, headers, lookup, inst)

    lookup = {
        Inst.BHROS: GeminiMapping,
        Inst.CIRPASS: Cirpass,
        Inst.F2: GeminiMapping,
        Inst.FLAMINGOS: Flamingos,
        Inst.ALOPEKE: Fox,
        Inst.ZORRO: Fox,
        Inst.GMOS: Gmos,
        Inst.GMOSS: Gmos,
        Inst.GMOSN: Gmos,
        Inst.GNIRS: GeminiMapping,
        Inst.GPI: GeminiMapping,
        Inst.GRACES: Graces,
        Inst.GSAOI: GeminiMapping,
        Inst.HOKUPAA: GeminiMapping,
        Inst.HRWFS: GeminiMapping,
        Inst.IGRINS: Igrins,
        Inst.MICHELLE: GeminiMapping,
        Inst.NICI: GeminiMapping,
        Inst.NIFS: GeminiMapping,
        Inst.NIRI: GeminiMapping,
        Inst.OSCIR: GeminiMapping,
        Inst.PHOENIX: GeminiMapping,
        Inst.TEXES: Texes,
        Inst.TRECS: GeminiMapping,
    }
    if inst in lookup:
        return lookup.get(inst)(storage_name, headers, metadata_lookup, inst)
    else:
        raise mc.CadcException(f'Mystery name {inst}.')
