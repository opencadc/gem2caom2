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
import traceback

from caom2 import Observation, CalibrationLevel, Chunk
from caom2 import TypedList, DerivedObservation
from caom2utils import WcsParser, update_artifact_meta
from caom2pipe import manage_composable as mc
from caom2pipe import caom_composable as cc
from caom2pipe import astro_composable as ac

import gem2caom2.external_metadata as em
import gem2caom2.obs_file_relationship as ofr
from gem2caom2.gem_name import GemName
from gem2caom2 import instruments, program_metadata
from gem2caom2.obs_metadata import json_lookup
from gem2caom2.util import Inst, COLLECTION, SCHEME


__all__ = ['GeminiMapping', 'APPLICATION']

APPLICATION = 'gem2caom2'


class GeminiMapping(cc.TelescopeMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def get_time_delta(self, ext):
        return em.current_instrument.get_time_delta(self._headers[ext])

    def get_calibration_level(self, ext):
        return em.current_instrument.get_calibration_level(
            self._storage_name.file_uri
        )

    def get_art_product_type(self, ext):
        return em.current_instrument.get_art_product_type(self._headers[ext])

    def get_data_product_type(self, ext):
        return em.current_instrument.get_data_product_type(self._headers[ext])

    def get_data_release(self, ext):
        """
        Determine the plane-level data release date.
        :return: The Plane release date, or None if not found.
        """
        # every instrument has a 'release' keyword in the JSON summary
        # not every instrument (Michelle) has a RELEASE keyword in
        # the appropriate headers
        result = json_lookup.get('release')
        if result is not None and result.startswith('0001'):
            # because obs id GN-2008A-Q-39-69-015
            result = result.replace('0001', '2001')
        return result

    def get_dec(self, ext):
        return em.current_instrument.get_dec(self._headers[ext])

    def get_exposure(self, ext):
        return em.current_instrument.get_exposure(self._headers[ext])

    def get_instrument_name(self, ext):
        return em.repair_instrument(json_lookup.get('instrument')).value

    def get_meta_release(self, ext):
        """
        Determine the metadata release date (Observation and Plane-level).
        :return: The Observation/Plane release date, or None if not found.
        """
        # make sure the metadata is for the correct plane/file
        # combination - this location happens to be the first function called
        # during blueprint evaluation, which is why reset is
        # called here
        file_id = ofr.remove_extensions(
            mc.CaomName(self._storage_name.file_uri).file_name
        )
        json_lookup.reset_index(file_id)

        # TODO - this is a bit different ROTFL
        if self._headers is None:
            # GenericParser, so no headers retrieved from archive.gemini.edu,
            # probably a 403 being returned by the site, assume proprietary
            meta_release = json_lookup.get('release')
        else:
            # DB 21-08-19
            # If PROP_MD is T, use JSON ‘release’ value for metadata release date.
            # If no PROP_MD present or value is F use the JSON ut_datetime value.
            prop_md = self._headers[ext].get('PROP_MD')
            if prop_md is None or prop_md is False or prop_md == 'F':
                meta_release = json_lookup.get('ut_datetime')
            else:
                meta_release = json_lookup.get('release')
        return meta_release

    def get_obs_intent(self, ext):
        return em.current_instrument.get_obs_intent(self._headers[ext])

    def get_obs_type(self, ext):
        return em.current_instrument.get_obs_type(self._headers[ext])

    def get_proposal_id(self, ext):
        """
        Determine the Proposal ID.

        :param header:  The FITS header for the current extension.
        :return: The proposal id from Gemini JSON metadata, or None if not found.
        """
        return json_lookup.get('program_id')

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
        return json_lookup.get('mode')

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
        cal_level = self.get_calibration_level(self._storage_name.file_uri)
        if cal_level in [
            CalibrationLevel.CALIBRATED,
            CalibrationLevel.PRODUCT,
            CalibrationLevel.ANALYSIS_PRODUCT,
        ]:
            comments = str(self._headers[ext].get('COMMENT'))
            result = fn(comments)
        return result

    def get_ra(self, ext):
        return em.current_instrument.get_ra(self._headers[ext])

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
        types = json_lookup.get('types')
        if 'NON_SIDEREAL' in types:
            return True
        else:
            return None

    def get_target_type(self, ext):
        return em.current_instrument.get_target_type()

    def get_time_function_val(self, ext):
        return em.current_instrument.get_time_function_val(self._headers[ext])

    def _get_data_label(self):
        return json_lookup.get('data_label')

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

        bp.set('Observation.type', 'get_obs_type(header)')
        bp.set('Observation.intent', 'get_obs_intent(header)')
        bp.set('Observation.metaRelease', 'get_meta_release(parameters)')
        bp.set('Observation.target.type', 'get_target_type(uri)')
        bp.set('Observation.target.moving', 'get_target_moving(header)')
        bp.set('Observation.proposal.id', 'get_proposal_id(header)')

        bp.clear('Observation.algorithm.name')
        bp.set('Observation.instrument.name', 'get_instrument_name(uri)')
        instrument = em.get_instrument(self._storage_name.file_uri)
        if instrument in [Inst.GMOSN, Inst.GMOSS, Inst.GMOS]:
            bp.set(
                'Observation.instrument.keywords',
                'get_provenance_keywords(uri)',
            )
        telescope = json_lookup.get('telescope')
        if telescope is not None or instrument is Inst.IGRINS:
            if telescope is not None and 'North' in telescope:
                x, y, z = ac.get_location(19.823806, -155.46906, 4213.0)
            else:
                x, y, z = ac.get_location(-30.240750, -70.736693, 2722.0)
            bp.set('Observation.telescope.geoLocationX', x)
            bp.set('Observation.telescope.geoLocationY', y)
            bp.set('Observation.telescope.geoLocationZ', z)

        bp.set('Plane.productID', self._storage_name.file_id)
        bp.set('Plane.dataProductType', 'get_data_product_type(header)')
        bp.set('Plane.calibrationLevel', 'get_calibration_level(uri)')
        bp.set('Plane.metaRelease', 'get_meta_release(parameters)')
        bp.set('Plane.dataRelease', 'get_data_release(header)')

        bp.set('Plane.provenance.name', 'Gemini Observatory Data')
        bp.set('Plane.provenance.project', 'Gemini Archive')
        # Add IMAGESWV for GRACES
        bp.add_fits_attribute('Plane.provenance.producer', 'IMAGESWV')
        bp.set_default('Plane.provenance.producer', 'Gemini Observatory')
        if instrument in [Inst.ALOPEKE, Inst.ZORRO]:
            bp.set(
                'Plane.provenance.reference',
                f'http://archive.gemini.edu/searchform/filepre='
                f'{self._storage_name.file_id}.fits',
            )
        elif instrument is not Inst.TEXES:
            data_label = self._get_data_label()
            bp.set(
                'Plane.provenance.reference',
                f'http://archive.gemini.edu/searchform/{data_label}',
            )

        if instrument is Inst.GRACES:
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

        bp.set('Artifact.productType', 'get_art_product_type(header)')
        bp.set('Artifact.contentChecksum', f'md5:{json_lookup.get("data_md5")}')
        bp.set('Artifact.contentLength', json_lookup.get('data_size'))
        bp.set('Artifact.contentType', 'application/fits')
        # always see the metadata, see the data only when it's public
        bp.set('Artifact.releaseType', 'data')
        bp.set('Artifact.uri', self._storage_name.file_uri)

        if instrument is Inst.CIRPASS:
            bp.set_default('Observation.telescope.name', 'Gemini-South')
        mode = json_lookup.get('mode')
        if not (
            instrument
            in [
                Inst.GPI,
                Inst.PHOENIX,
                Inst.HOKUPAA,
                Inst.OSCIR,
                Inst.BHROS,
                Inst.TRECS,
            ]
            or (
                instrument is Inst.GRACES
                and mode is not None
                and mode != 'imaging'
            )
        ):
            bp.configure_position_axes((1, 2))

        if instrument is Inst.FLAMINGOS:
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
        if instrument in [Inst.ALOPEKE, Inst.ZORRO]:
            bp.clear('Chunk.time.axis.function.naxis')
            bp.add_fits_attribute('Chunk.time.axis.function.naxis', 'NAXIS3')
            bp.set_default('Chunk.time.axis.function.naxis', 1)

        bp.set('Chunk.time.axis.function.delta', 'get_time_delta(header)')
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

        instrument = Inst(observation.instrument.name)

        # processed files
        if cc.is_composite(self._headers) and not isinstance(
            observation, DerivedObservation
        ):
            observation = self._update_composite(observation, instrument)

        if instrument in [Inst.MICHELLE, Inst.GNIRS]:
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
                    f'{instrument}: no image data for '
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
                    json_lookup.reset_index(file_id)
                    processed = ofr.is_processed(caom_name.file_name)
                    if instrument in [
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
                            artifact, instrument, observation.observation_id
                        )

                    delete_these_parts = []
                    for part in artifact.parts:

                        if part == '2' and instrument is Inst.GPI:
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
                            x = instruments.instrument_factory(instrument)
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
                            mode = json_lookup.get('mode')
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
                        or instrument is Inst.TEXES
                    ) and 'jpg' not in caom_name.file_name:
                        # not the preview artifact
                        if plane.provenance is not None:
                            if instrument is not Inst.GRACES:
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
                f'{instrument}'
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

    def _update_position_from_zeroth_header(self, artifact, instrument, obs_id):
        """Make the 0th header spatial WCS the WCS for all the
        chunks."""
        primary_header = self._headers[0]

        # naxis values are only available from extensions

        primary_header['NAXIS1'] = self._headers[-1].get('NAXIS1')
        primary_header['NAXIS2'] = self._headers[-1].get('NAXIS2')

        primary_chunk = Chunk()
        types = json_lookup.get('types')
        ra = self.get_ra(0)
        ra_json = json_lookup.get('ra')
        dec_json = json_lookup.get('dec')
        if ('AZEL_TARGET' in types and ra is None) or (
            instrument is Inst.GNIRS and ra_json is None and dec_json is None
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

            self._logger.info(f'{instrument}: Spatial WCS is None for {obs_id}')
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

    def _update_composite(self, obs, instrument):
        result = None
        if instrument is Inst.TRECS:
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
