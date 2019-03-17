# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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

from caom2pipe import manage_composable as mc
from caom2pipe import execute_composable as ec
from gem2caom2.gem_name import GemName


class GeminiObsMetadata(object):
    """A place to hold access to output from multiple jsonsummary
    queries.

    Hold the query results for all files associated with an observation.
    Use the 'add' method to add a single jsonsummary query result.
    Use the 'reset_index' method to have the 'get' method look up the
    results associated with a particular file_id.
    """

    def __init__(self):
        # a dictionary of all the jsonsummary results
        self.lookup = {}
        # which dictionary entry is of current lookup interest
        self.current = None
        # the json summary results are a list, track which entry in the
        # list has the information for a particular file_id
        self.index = -1

    def add(self, metadata, file_id):
        self.lookup[file_id] = metadata
        self._reset_index(file_id)

    def get(self, lookup):
        return mc.response_lookup(self.current[self.index], lookup)

    def reset_index(self, uri):
        file_id = GemName.remove_extensions(ec.CaomName(uri).file_name)
        self._reset_index(file_id)

    def _reset_index(self, file_id):
        if file_id not in self.lookup:
            raise mc.CadcException(
                'ObsMetadata: Mystery file id {}'.format(file_id))
        self.current = self.lookup[file_id]
        self.index = self._get_index(file_id)

    def _get_index(self, file_id):
        result = -1
        for index, value in enumerate(self.current):
            indexed_f_name = mc.response_lookup(value, 'filename')
            if indexed_f_name is not None:
                temp = GemName.remove_extensions(indexed_f_name)
                if temp == file_id:
                    result = index
                    break
        if result == -1:
            # TODO - set obs id?
            raise mc.CadcException(
                'JSON Summary: unrecognized file_id {} in obs_id {}'.format(
                    file_id, ''))
        return result

# DB -  12-03-19 - Processed Observation Guidance

# GMOS has some processed ‘fringe’ datasets.  One for Gemini South has file
# name rgS20161227S0051_fringe.fits, with datalabel GS-CAL20161227-5-001.
# But note that same datalabel is used for the raw image S20161227S0051.fits.
# But since the ‘fringe’ version is a composite observation (with NCOMBINE
# and IMCMB### keywords) the processed observation should NOT overwrite the
# raw simple observation and isn’t a second plane.  Not sure how you handle
# this type of thing.  Could add ‘-rg-fringe’ to the datalabel.
#
# But here’s another one for Gemini South:  file S20131007S0067_fringe is
# also a coadded fringe frame with datalabel.  GS-CAL20131007-900-067.
#
# And a Gemini North fringe file:  filename = N20120320S0328_stack_fringe;
# datalabel = GN-CAL20120320-900-328.  This one has NCOMBINE=5 but does NOT
# have IMCMB### keywords so will have to be stored as a simple observation.
#
# And another Gemini North file with different naming convention:
# filename = N20110927S0170_fringe; datalabel = GN-CAL20110927-900-170.
#
# Here’s one in Paul’s list,
# rgS20180122S0236_fringe = GS-CAL20180122-1-001, but when I search for
# that filename in the Gemini archive interface the datalabel shown is
# GS-CAL20180122-1-001-rg-fringe.
#
# And finally, mrgS20181016S0184_fringe / GS-CAL20181016-5-001.
# The prefix ‘m’ means ‘GMOSAIC’ has been run - all of the image extensions
# have been combined into a single large image/extension.
#
# And GMOS processed science images (likely used to make mask files).
# Samples:
#
# mrgS20160901S0122_add = GS-2016B-Q-72-23-001-MRG-ADD
# mrgN20160311S0691_add = GN-2016A-Q-68-46-001-MRG-ADD
# mfrgS20160310S0154_add = GS-2016A-Q-7-175-001-MFRG-ADD
#
# Ideally, since these are co-added science observations the temporal WCS
# should be sampled to include the time info for the individual subexposures
# if identified.

# Some datalabels are duplicated, some datalabels do not exist, and at
# least the one example above appears to have a different datalabel from
# what Paul provided us in his file.

# TReCS processed data are much more straightforward (I think GMOS data make
# up well over 50% of all processed datasets, and co-added biases are about
# 80% of the processed GMOS datasets).  There are about 9,500 processed
# TReCS datasets.  ALL should be added as another plane with cal. level = 2.
# Note json ‘reduction’ = ‘prepared’ for these files apparently.
#
# filename:  rS20121030S0136   datalabel:  GS-2012B-Q-90-366-003-R
#
# So the one above would be added as a CAL=2 plane for the observation
# GS-2012B-Q-90-366-003.  Product ID could be the new datalabel.

# Another type of processed GMOS science image:
#
# rgS20100316S0366 = GS-2010A-Q-36-6-358-RG
#
# I think these are relatively rare cases of GMOS-N/S images being obtained
# with only one CCD so no need for “mosaic”ing the 3 CCD’s into one image.
# No ‘add’ in the name/label so these are additional planes to the raw
# observation.

# NIRI processed files (conventions sometimes match GMOS-N/S) are composite
# flat and dark observations.  json ‘reduction’ = ‘prepared’ for these as
# well (also the case for GMOS-N/S):
#
# N20140313S0072_flat = GN-2013B-Q-75-163-011_STACK
# N20150804S0348_dark GN-2015B-Q-53-138-061_STACK
# or
# N20070819S0339_dark GN-2007B-Q-107-150-004_DARK
# N20130404S0512_flat GN-2013A-Q-63-54-051_FLAT
#
# They don’t always have IMCMB### keywords to identify members. (edited)
# It looks like F2 also uses the same naming conventions for flat/dark
# exposures as NIRI:
#
# S20140124S0039_dark = GS-2013B-Q-16-277-019_STACK
# or
# S20141129S0331_dark = GS-CAL20141129-1-001_DARK

# As far as I can tell all of NIRI, F2, GNIRS and GSAOI use the same
# flat/STACK, dark/STACK, flat/FLAT and dark/DARK filename/datalabel
# naming conventions.  There are NOT very many of these processed datasets
# in the archive.

# PHOENIX.  Basically anything that has a filename that starts with “P|p”
# is a processed dataset.  Little header information about processing,
# suffix of file name appears to be freeform.  I think the only thing we
# can do is add these as planes to whatever observation is identified in
# their datalabel except possibly for some files with “COMB” or “FLAT” in
# the suffix.  These might have NCOMBINE and IMCMB### headers to identify
# members in which case they could be made composite observations.
# Filename examples:
#
# p2004may20_0048_FLAT = GS-CAL20040520-7-0048-FLAT-P
# p2004may19_0255_COMB = GS-2004A-Q-6-27-0255-COMB-P
# P2003JAN14_0148_DARK = GS-CAL20030114-7-0148
# P2002FEB03_0045_DARK10SEC = GS-CAL20020203-4-0045
# P2002DEC02_0161_SUB = GS-2002B-Q-22-13-0161
# P2002DEC02_0075_SUB.0001 = GS-CAL20021202-3-0075

# And I just noticed that for some reason all PHOENIX files obtained in 2016
# May have a filename prefix ‘c’ with other text identifying the obstype.
# e.g. c2016may18_sci061.fits, c2016may18_acq047.fits.  These are not
# processed files.  Only happens in 2016 May - PHOENIX is a visiting
# instrument that reappears at Gemini-South every so often.

# DB - 04-03-19  TODO
# Example of GMOS composites:  any datasets with datalabels like
# GS-CAL20190301-4-046-g-bias, GS-CAL20181219-4-021-g-flat and
# GS-CAL20130103-3-001-rg-fringe are composite processed bias/flat/fringe
# images.  Corresponding file names are gS20190301S0556_bias.fits,
# gS20181219S0216_flat.fits, and rgS20130103S0098_FRINGE.fits.
# (With some naming inconsistencies perhaps.)   It also won’t always be
# possible to determine the members since the information isn’t always in
# the headers.
# GMOS:  I don’t think there are any GMOS processed science observations
# that are ‘composites’ but there will be some images with one or two
# processed planes.  e.g. datalable GS-2010A-Q-36-5-246-rg, file name
# rgS20100212S0301.fits.  Basically any Gemini image with a ‘-[a-z|A-Z]’
# suffix to the datalabel is a processed dataset.
#
# File names can have a suffix AND a prefix.

