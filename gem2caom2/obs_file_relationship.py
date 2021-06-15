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

"""
DB -  12-03-19 - Processed Observation Guidance

GMOS has some processed ‘fringe’ datasets.  One for Gemini South has file
name rgS20161227S0051_fringe.fits, with datalabel GS-CAL20161227-5-001.
But note that same datalabel is used for the raw image S20161227S0051.fits.
But since the ‘fringe’ version is a composite observation (with NCOMBINE
and IMCMB### keywords) the processed observation should NOT overwrite the
raw simple observation and isn’t a second plane.  Not sure how you handle
this type of thing.  Could add ‘-rg-fringe’ to the datalabel.

But here’s another one for Gemini South:  file S20131007S0067_fringe is
also a coadded fringe frame with datalabel.  GS-CAL20131007-900-067.

And a Gemini North fringe file:  filename = N20120320S0328_stack_fringe;
datalabel = GN-CAL20120320-900-328.  This one has NCOMBINE=5 but does NOT
have IMCMB### keywords so will have to be stored as a simple observation.

And another Gemini North file with different naming convention:
filename = N20110927S0170_fringe; datalabel = GN-CAL20110927-900-170.

Here’s one in Paul’s list,
rgS20180122S0236_fringe = GS-CAL20180122-1-001, but when I search for
that filename in the Gemini archive interface the datalabel shown is
GS-CAL20180122-1-001-rg-fringe.

And finally, mrgS20181016S0184_fringe / GS-CAL20181016-5-001.
The prefix ‘m’ means ‘GMOSAIC’ has been run - all of the image extensions
have been combined into a single large image/extension.

And GMOS processed science images (likely used to make mask files).
Samples:

mrgS20160901S0122_add = GS-2016B-Q-72-23-001-MRG-ADD
mrgN20160311S0691_add = GN-2016A-Q-68-46-001-MRG-ADD
mfrgS20160310S0154_add = GS-2016A-Q-7-175-001-MFRG-ADD

Ideally, since these are co-added science observations the temporal WCS
should be sampled to include the time info for the individual subexposures
if identified.

Some datalabels are duplicated, some datalabels do not exist, and at
least the one example above appears to have a different datalabel from
what Paul provided us in his file.

TReCS processed data are much more straightforward (I think GMOS data make
up well over 50% of all processed datasets, and co-added biases are about
80% of the processed GMOS datasets).  There are about 9,500 processed
TReCS datasets.  ALL should be added as another plane with cal. level = 2.
Note json ‘reduction’ = ‘prepared’ for these files apparently.

filename:  rS20121030S0136   datalabel:  GS-2012B-Q-90-366-003-R

So the one above would be added as a CAL=2 plane for the observation
GS-2012B-Q-90-366-003.  Product ID could be the new datalabel.

Another type of processed GMOS science image:

rgS20100316S0366 = GS-2010A-Q-36-6-358-RG

I think these are relatively rare cases of GMOS-N/S images being obtained
with only one CCD so no need for “mosaic”ing the 3 CCD’s into one image.
No ‘add’ in the name/label so these are additional planes to the raw
observation.

NIRI processed files (conventions sometimes match GMOS-N/S) are composite
flat and dark observations.  json ‘reduction’ = ‘prepared’ for these as
well (also the case for GMOS-N/S):

N20140313S0072_flat = GN-2013B-Q-75-163-011_STACK
N20150804S0348_dark GN-2015B-Q-53-138-061_STACK
or
N20070819S0339_dark GN-2007B-Q-107-150-004_DARK
N20130404S0512_flat GN-2013A-Q-63-54-051_FLAT

They don’t always have IMCMB### keywords to identify members. (edited)
It looks like F2 also uses the same naming conventions for flat/dark
exposures as NIRI:

S20140124S0039_dark = GS-2013B-Q-16-277-019_STACK
or
S20141129S0331_dark = GS-CAL20141129-1-001_DARK

As far as I can tell all of NIRI, F2, GNIRS and GSAOI use the same
flat/STACK, dark/STACK, flat/FLAT and dark/DARK filename/datalabel
naming conventions.  There are NOT very many of these processed datasets
in the archive.

PHOENIX.  Basically anything that has a filename that starts with “P|p”
is a processed dataset.  Little header information about processing,
suffix of file name appears to be freeform.  I think the only thing we
can do is add these as planes to whatever observation is identified in
their datalabel except possibly for some files with “COMB” or “FLAT” in
the suffix.  These might have NCOMBINE and IMCMB### headers to identify
members in which case they could be made composite observations.
Filename examples:

p2004may20_0048_FLAT = GS-CAL20040520-7-0048-FLAT-P
p2004may19_0255_COMB = GS-2004A-Q-6-27-0255-COMB-P
P2003JAN14_0148_DARK = GS-CAL20030114-7-0148
P2002FEB03_0045_DARK10SEC = GS-CAL20020203-4-0045
P2002DEC02_0161_SUB = GS-2002B-Q-22-13-0161
P2002DEC02_0075_SUB.0001 = GS-CAL20021202-3-0075

And I just noticed that for some reason all PHOENIX files obtained in 2016
May have a filename prefix ‘c’ with other text identifying the obstype.
e.g. c2016may18_sci061.fits, c2016may18_acq047.fits.  These are not
processed files.  Only happens in 2016 May - PHOENIX is a visiting
instrument that reappears at Gemini-South every so often.

DB - 04-03-19  TODO
Example of GMOS composites:  any datasets with datalabels like
GS-CAL20190301-4-046-g-bias, GS-CAL20181219-4-021-g-flat and
GS-CAL20130103-3-001-rg-fringe are composite processed bias/flat/fringe
images.  Corresponding file names are gS20190301S0556_bias.fits,
gS20181219S0216_flat.fits, and rgS20130103S0098_FRINGE.fits.
(With some naming inconsistencies perhaps.)   It also won’t always be
possible to determine the members since the information isn’t always in
the headers.
GMOS:  I don’t think there are any GMOS processed science observations
that are ‘composites’ but there will be some images with one or two
processed planes.  e.g. datalabel GS-2010A-Q-36-5-246-rg, file name
rgS20100212S0301.fits.  Basically any Gemini image with a ‘-[a-z|A-Z]’
suffix to the datalabel is a processed dataset.

File names can have a suffix AND a prefix.

DB - 18-03-19
The processed PHOENIX file GS-2004A-Q-6-27-0255-COMB-P shows a problem
with PHOENIX filters that I wasn’t aware of.  The actual filter in the
beam can be identified with either the value of the FIL_POS keyword or
the CVF_POS.  In this case FIL_POS contains the string ‘open’ and should
be ignored and the CVF_POS value used to look up energy info.  Maybe add
what is likely the first file added in this ‘comb’ file to the test:
GS-2004A-Q-6-27-0255.

P2002DEC02_0161_SUB.fits has a ‘regular’ datalabel,
‘GS-2002B-Q-22-13-0161’.   Can you add the latter unprocessed file as a
test to see if the ‘sub’ file is added as another plane? (The file name
is 2002dec02_0161.fits)

GMOS:  could you add GS-CAL20181016-5-001 as a test unprocessed file?
Then mrgS20181016S0184_fringe.fits should become a second processed
plane.   I guess I should have given you the unprocessed versions of all
of the non-co-added examples.

<file id>_flat_pasted is a different observation than <file id>_flat.


DB - 02-06-20

Most TReCS files are being processed correctly except for the 'composite'
algorithm. The exceptions might be those with data labels ending in -G. I may
have only found examples with -R suffixes previously. The latter are correctly
combined with the unprocessed observations (again, except for 'composite')
whereas the -G versions show up as distinct observations.
"""

import collections
import logging
import re

from datetime import datetime
from datetime import timedelta

from caom2pipe import manage_composable as mc

from gem2caom2 import gem_name


__all__ = ['GemObsFileRelationship', 'FILE_NAME', 'repair_data_label']

FILE_NAME = '/app/data/from_paul.txt'


class GemObsFileRelationship(object):
    """A class to hold and access the content of the observation ID/file id
    information that is provided to CADC from Gemini, while also adhering to
    the file-to-observation guidance laid out in the module comment.

    It's made into a class, because the information is useful from both
    the gem2caom2 repo, for identifying provenance relationships, and from
    the gemHarvester2Caom2, for supporting list_observations queries.
    """

    def __init__(self):

        # id_list structure: a dict, keys are Gemini observation IDs, values
        # are a set of associated file names. This structure supports the
        # get_observation query.

        self.id_list = collections.defaultdict(list)

        # time_list structure: a dict, keys are last modified time,
        # values are a set of observation IDs as specified from Gemini
        # with that last modified time. This structure supports the
        # time-bounded queries of the Harvester.

        self.time_list = {}

        # name_list structure: a dict, keys are file ids, values are Gemini
        # observation IDs. This structure supports queries by gem2caom2
        # for determining provenance information for planes and
        # observations.

        self.name_list = collections.defaultdict(list)

        self.logger = logging.getLogger(__name__)
        self._initialize_content(FILE_NAME)

    def _initialize_content(self, fqn):
        """Initialize the internal data structures that represents the
        query list from the Gemini Science Archive.
        """
        result = self._read_file(fqn)
        # result row structure:
        # 0 = data label
        # 1 = timestamp
        # 2 = file name
        temp_content = {}
        logging.info('Progress - file read ....')
        for ii in result:
            # re-organize to be able to answer list_observations queries
            ol_key = mc.make_seconds(ii[1])
            if ol_key in temp_content:
                if ii[0] not in temp_content[ol_key]:
                    temp_content[ol_key].append(ii[0])
            else:
                temp_content[ol_key] = [ii[0]]
            # re-organize to be able to answer get_observation queries
            self.id_list[ii[0]].append(ii[2])
            file_id = gem_name.GemName.remove_extensions(ii[2])
            self.name_list[file_id].append([ii[0], ol_key])

        # this structure means an observation ID occurs more than once with
        # different last modified times
        self.time_list = collections.OrderedDict(
            sorted(temp_content.items(), key=lambda t: t[0])
        )
        self.logger.info('Observation list initialized in memory.')

    def _read_file(self, fqn):
        """Read the .txt file from Gemini, and make it prettier ...
        where prettier means stripping whitespace, query output text, and
        making an ISO 8601 timestamp from something that looks like this:
        ' 2018-12-17 18:19:27.334144+00 '

        or this:
        ' 2018-12-17 18:19:27+00 '

        :return a list of lists, where the inner list consists of an
            observation ID, a last modified date/time, and a file name.

        File structure indexes:
        0 == data label
        1 == file name
        3 == last modified date/time
        """
        results = []
        try:
            with open(fqn) as f:
                for row in f:
                    temp = row.split('|')
                    if len(temp) > 1 and 'data_label' not in row:
                        time_string = temp[3].strip().replace(' ', 'T')
                        if '/' in temp[0]:
                            if 'MBIAS' in temp[0]:
                                temp[0] = temp[0].replace('BIAS/MBIAS/', '')
                            elif 'PETRO' in temp[0]:
                                temp[0] = temp[0].replace(
                                    '-/NET/PETROHUE/DATAFLOW/', ''
                                )
                            elif '12CD' in temp[0] or 'EXPORT/HOME' in temp[0]:
                                temp[0] = temp[0].split('/', 1)[0]
                            else:
                                logging.warning(
                                    f'Mystery data label {temp[0]}'
                                )
                        elif '?' in temp[0]:
                            if 'GS-2002A-DD-1-?' in temp[0]:
                                temp[0] = temp[0].replace('?', '11')
                            else:
                                logging.warning(
                                    f'Mystery data label {temp[0]}'
                                )
                        elif '"' in temp[0]:
                            temp[0] = temp[0].replace('"', '')
                        if len(temp[0].strip()) > 1:
                            results.append(
                                [temp[0].strip(), time_string, temp[1].strip()]
                            )
                        else:
                            # no data label in the file, so use the file name
                            results.append(
                                [temp[1].strip(), time_string, temp[1].strip()]
                            )

        except Exception as e:
            self.logger.error(f'Could not read from csv file {fqn}')
            raise mc.CadcException(e)
        return results

    def subset(self, start=None, end=None, maxrec=None):
        if start is not None and end is not None:
            temp = self._subset(start.timestamp(), end.timestamp())
        elif start is not None:
            temp = self._subset(start.timestamp(), datetime.now().timestamp())
        elif end is not None:
            temp = self._subset(0, end.timestamp())
        else:
            temp = self._subset(0, datetime.now().timestamp())
        if maxrec is not None:
            temp = temp[:maxrec]
        return temp

    def _subset(self, start_s, end_s):
        """Get only part of the observation list, limited by timestamps."""
        self.logger.debug(
            f'Timestamp endpoints are between {start_s} and {end_s}.'
        )
        temp = []
        for ii in self.time_list:
            if start_s <= ii <= end_s:
                for jj in self.time_list[ii]:
                    dt = datetime.fromtimestamp(ii).isoformat(
                        timespec="milliseconds"
                    )
                    temp.append(f'{gem_name.COLLECTION} {jj} {dt}')
            if ii > end_s:
                break
        return temp

    def get_file_names(self, obs_id):
        """Given an obs id, return the list of file names that make up
        the observation."""
        if obs_id in self.id_list:
            temp_str = ''.join(self.id_list[obs_id])
            if re.search(r'_BIAS|_FLAT', temp_str) is not None:
                # all of this is to handle the approximately 5000 cases where
                # the same data label refers to different file names, where
                # the difference is in the case of the file name only
                temp_set = set([ii.upper() for ii in self.id_list[obs_id]])
                if len(temp_set) != len(self.id_list[obs_id]):
                    temp_list = []
                    for f_name in self.id_list[obs_id]:
                        x = self._check_duplicate(f_name.replace('.fits', ''))
                        temp_list.append(f'{x}.fits')
                    return list(set(temp_list))
            else:
                return self.id_list[obs_id]
        else:
            return None

    def get_obs_id(self, file_id):
        checked = self._check_duplicate(file_id)
        if checked in self.name_list:
            # structure of the entry is ['obs id', timestamp], so return
            # only the obs_id of the first entry
            return self.name_list[checked][0][0]
        else:
            return None

    def get_timestamp(self, file_id):
        checked = self._check_duplicate(file_id)
        if checked in self.name_list:
            temp = self.name_list[checked]
            return temp[0][1]
        else:
            return timedelta()

    def get_max_timestamp(self):
        return list(self.time_list.keys())[-1]

    def repair_data_label(self, file_id):
        """For processed files, try to provide a consistent naming pattern,
        because data labels aren't unique within Gemini, although the files
        they refer to are, and can be in different CAOM Observations.

        Take the prefixes and suffixes on the files, that indicate the type of
        processing, and append them in upper case, to the data label, for
        uniqueness.

        DB - 07-03-19
        TEXES Spectroscopy

        Some special code will be needed for datalabels/planes.  There are no
        datalabels in the FITS header.  json metadata (limited) must be
        obtained with URL like
        https://archive.gemini.edu/jsonsummary/canonical/filepre=TX20170321_flt.2507.fits.
        Use TX20170321_flt.2507 as datalabel.  But NOTE:  *raw.2507.fits and
        *red.2507.fits are two planes of the same observation. I’d suggest we
        use ‘*raw*’ as the datalabel and ‘*red*’ or ‘*raw*’ as the appropriate
        product ID’s for the science observations.  The ‘flt’ observations do
        not have a ‘red’ plane.  The json document contains ‘filename’ if
        that’s helpful at all.  The ‘red’ files do not exist for all ‘raw’
        files.
        """
        if file_id in self.name_list:
            repaired = self.name_list[file_id][0][0]
            # if the data label is missing, the file name, including
            # extensions, is treated as the data label, so get rid of .fits
            repaired = gem_name.GemName.remove_extensions(repaired)
            repaired = repair_data_label(file_id, repaired)
        else:
            logging.warning(
                f'File name {file_id} not found in the Gemini list.'
            )
            repaired = file_id
        return repaired

    def _check_duplicate(self, file_id):
        """There are data labels in the Gemini-supplied file, where the
        only difference in the related file name is the case of the 'bias'
        or 'flat' text. So look for that as well, when checking for
        membership.

        There are approximately 5500 entries of this duplicate sort - too
        many to fix manually, but also too many to ignore.
        """
        if re.search(r'_BIAS|_FLAT', file_id) is not None:
            duplicate_check = re.sub(
                '_FLAT',
                '_flat',
                re.sub('_BIAS', '_bias', file_id),
            )
            if duplicate_check in self.name_list:
                logging.warning(f'Replacing {file_id} with {duplicate_check}')
                file_id = duplicate_check
        return file_id


def get_prefix(file_id):
    if file_id.startswith(('p', 'P')):
        if '_FLAT' in file_id or '_COMB' in file_id:
            prefix = 'P'
        else:
            prefix = ''
    elif file_id.startswith(('TX', 'ag', 'c')):
        prefix = ''
    elif -1 < file_id.find('GN') < 14:
        prefix = file_id.split('GN', 1)[0]
    elif -1 < file_id.find('N') < 14:
        prefix = file_id.split('N', 1)[0]
    elif 'GS' in file_id:
        prefix = file_id.split('GS', 1)[0]
    elif 'S' in file_id:
        prefix = file_id.split('S', 1)[0]
    else:
        logging.warning(f'Unrecognized file_id pattern {file_id}')
        prefix = ''
    return prefix


def get_suffix(file_id, data_label):
    temp = []
    suffix = []
    if '-' in file_id:
        temp = file_id.split('-')[:1]
    elif '_' in file_id:
        if file_id.startswith(('p', 'P')):
            if (
                '_FLAT' in file_id
                or '_COMB' in file_id
                or '_flat' in file_id
                or '_comb' in file_id
            ):
                temp = file_id.split('_')[2:]
        elif file_id.startswith('TX2'):
            # when the data label is the file id, fix every
            # data label except flats
            if not data_label.startswith(file_id):
                if '_flt' in file_id.lower():
                    temp = ['flt']
        else:
            temp = file_id.split('_')[1:]
    if data_label.endswith('-G') and (
        file_id.startswith('rS') or file_id.startswith('rN')
    ):
        # DB 16-06-20
        # I think the ‘g’ prefix is used a little inconsistently.  It is
        # supposed to be set whenever the IRAF GPREPARE is executed and I
        # think this is normally done at the start of ALL Gemini processing.
        # The ‘r’ prefix should appear in all processed files.
        temp.append('g')
    for ii in temp:
        if re.match('[a-zA-Z]+', ii) is not None:
            suffix.append(ii)
    return suffix


def get_removals(file_id, repaired):
    removals = []
    if file_id.startswith('TX2'):
        # when the data label is the file id, fix every
        # data label except flats
        if repaired.startswith(file_id):
            if '_flt' not in file_id.lower():
                removals = ['_raw', '_red', '_sum']
        else:
            if '_flt' not in file_id.lower():
                removals = ['raw', 'red', 'sum']
    return removals


def is_processed(file_name):
    """Try to determine if a Gemini file is processed, based on naming
    patterns."""
    result = True
    file_id = gem_name.GemName.remove_extensions(file_name)
    # ALOPEKE file id ends with 'r' or 'b', so avoid checking that letter
    if file_id.startswith(('S', 'N', 'GN', 'GS', 'c', 'abu')):
        if file_id.endswith(
            ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
        ):
            result = False
        if file_id[:15].endswith(('b', 'r')):
            result = False
    elif file_id.startswith(('2', '02', '01')):
        result = False
    # TEXES naming patterns
    elif file_id.startswith('TX2') and '_raw' in file_id:
        result = False
    # OSCIR file naming pattern
    elif (
        file_id.startswith('r')
        and re.match('r\\w{7}_\\d{3}', file_id, flags=re.ASCII) is not None
    ):
        result = False
    return result


def repair_data_label(file_name, data_label):
    """For processed files, try to provide a consistent naming pattern,
    because data labels aren't unique within Gemini, although the files
    they refer to are, and can be in different CAOM Observations.

    Take the prefixes and suffixes on the files, that indicate the type of
    processing, and append them in upper case, to the data label, for
    uniqueness.

    DB - 07-03-19
    TEXES Spectroscopy

    Some special code will be needed for datalabels/planes.  There are no
    datalabels in the FITS header.  json metadata (limited) must be
    obtained with URL like
    'https://archive.gemini.edu/jsonsummary/canonical/filepre=
     TX20170321_flt.2507.fits.'

    Use TX20170321_flt.2507 as datalabel.  But NOTE:  *raw.2507.fits and
    *red.2507.fits are two planes of the same observation. I’d suggest we
    use ‘*raw*’ as the datalabel and ‘*red*’ or ‘*raw*’ as the appropriate
    product ID’s for the science observations.  The ‘flt’ observations do
    not have a ‘red’ plane.  The json document contains ‘filename’ if
    that’s helpful at all.  The ‘red’ files do not exist for all ‘raw’
    files.
    """
    # if the data label is missing, the file name, including
    # extensions, is treated as the data label, so get rid of .fits
    file_id = gem_name.GemName.remove_extensions(file_name)
    repaired = data_label
    if is_processed(file_id) or file_id.startswith('TX2'):
        if not file_id.startswith('TX2'):
            repaired = repaired.split('_')[0]

        prefix = get_prefix(file_id)
        suffix = get_suffix(file_id, repaired)
        removals = get_removals(file_id, repaired)
        if len(prefix) > 0:
            removals = [prefix] + suffix
        else:
            removals = removals + suffix

        for ii in removals:
            # rreplace
            temp = repaired.rsplit(ii, 1)
            repaired = ''.join(temp)
            temp = repaired.rsplit(ii.upper(), 1)
            repaired = ''.join(temp)
            repaired = repaired.rstrip('-')
            repaired = repaired.rstrip('_')

        # DB - 18-03-19
        # Basically any ‘mfrg’, ‘mrg’, or ‘rg’ file WITHOUT ‘add’
        # in the datalabel or name is a processed version of a raw
        # file without the datalabel suffix (filename prefix)
        #
        # DB - 19-07-31
        # Assuming all such processed ‘arc’ files are processed
        # identically then these should be different planes in the
        # same observation.  The “_arc” file has had some basic
        # processing carried out.
        #
        # r<file name> should be another plane of the same
        # observation.
        #
        # DB - 02-06-20
        # r<file name> should NOT be composites. I believe processed TReCS
        # files in Gemini's archive are derived by combining the
        # NNODSETS x NSAVSETS contained within a single unprocessed image
        # into a simpler image array.
        #
        # Most TReCS files are being processed correctly except for the
        # 'composite' algorithm. The exceptions might be those with data
        # labels ending in -G. The -G versions should show up as an additional
        # plane in a single observation.
        #
        # SGo - this means make the data labels the same
        if (
            (
                ('mfrg' == prefix or 'mrg' == prefix or 'rg' == prefix)
                and (
                    not (
                        'add' in suffix
                        or 'ADD' in suffix
                        or 'fringe' in suffix
                        or 'FRINGE' in suffix
                    )
                )
            )
            or ('arc' in suffix or 'ARC' in suffix)
            or ('r' == prefix or 'R' == prefix)
        ):
            prefix = ''
            suffix = []

        if (
            prefix == ''
            and len(suffix) == 1
            and ('FRINGE' in suffix or 'fringe' in suffix)
        ):
            suffix = []

        if len(prefix) > 0:
            if f'-{prefix.upper()}' not in repaired:
                repaired = f'{repaired}-{prefix.upper()}'

        for ii in suffix:
            if f'-{ii.upper()}' not in repaired:
                repaired = f'{repaired}-{ii.upper()}'
    else:
        repaired = file_id if repaired is None else repaired
    return repaired
