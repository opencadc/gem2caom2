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

import logging
import re


__all__ = ['repair_data_label', 'remove_extensions']


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
    file_id = remove_extensions(file_name)
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

    ALOPEKE/ZORRO::

    DB 31-08-20
    DATALAB can NOT be used for the CAOM2 Observation ID since it appears that
    the DATALAB value is identical for all files obtained for a single
    program. e.g. if the program ID is GN-2020A-DD-115 then the DATALAB value
    is always GN-2020A-DD-115-0-0.

    Instead, use the root of the filename as the observation ID.  e.g.
    N20200819A0003r.fits and N20200819A0003b.fits are two files generated from
    a single observation (r = red channel, b = blue channel).  Use
    N20200819A0003 as the observation ID with two planes given by the two
    colours of data.

    DB 01-09-20
    Gemini has kludged the headers so that every observation for a single
    program has the same DATALAB in the header.  This is what we usually use
    for the observation ID.  Each single ‘observation’ actually produces two
    files (not a single MEF file) for the red and blue channels so to me it
    would make the most sense to group these two files as a single observation
    with two artifacts given by uri’s pointing to the two files.  And this is
    a single plane, correct?

    PD 01-09-20
    What is the meaning of red and blue channels? different energy bands?

    DB 02-09-20
    Yes.  there’s a dichroic that directs light shortward of 675nm to one
    detector (through one of several possible filters) and light longward of
    675nm to a second detector (through another filter).   But instead of
    generating a single MEF file they generate two files, e.g.
    N20191219A0004b.fits and N20191219A0004r.fits.

    PD 02-09-20
    This seems very much like MACHO... if those two files are images in the
    normal sense then it could make sense to create separate planes with
    dataProductType = image that end up with the correct (distinct) energy
    metadata. It is OK for an observation to create two sibling products and
    two planes probably captures the goal of this instrument/observing mode
    more directly.

    DB 21-07-21
    IGRINS modify DATALAB values to give (hopefully) unique observation IDs
    - use an observation ID of GS-2020B-Q-315-23-1104 instead of
    GS-2020B-Q-315-23-0 for file SDCH_20201104_0023.fits, by grabbing the
    MMDD from the file name (1104 in this case) and replacing the trailing
    -0 with -MMDD.
    """
    # if the data label is missing, the file name, including
    # extensions, is treated as the data label, so get rid of .fits
    logging.debug(
        f'Begin repair_data_label with file {file_name} and data label '
        f'{data_label}.'
    )
    file_id = remove_extensions(file_name)
    repaired = data_label if data_label else ''
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
    elif file_id.endswith('r') or file_id.endswith('b'):
        # Alopeke/Zorro files, data_label is the file_id minus the
        # channel indicator
        repaired = file_id[:-1]
    elif file_id.startswith('SDC'):
        # IGRINS
        file_id_bits = file_id.split('_')
        data_label_good_bits = data_label.rsplit('-0', 1)
        repaired = f'{data_label_good_bits[0]}-{file_id_bits[1][4:]}'
    else:
        repaired = file_id if repaired is None else repaired
    logging.debug(
        f'End repair_data_label with file {file_name} and data label '
        f'{repaired}.'
    )
    return repaired


def remove_extensions(name):
    """How to get the file_id from a file_name."""
    # Note the .gz extension is on some TRECS files, not that it is
    # an accepted GEMINI extension
    return (
        name.replace('.fits', '')
            .replace('.bz2', '')
            .replace('.header', '')
            .replace('.jpg', '')
            .replace('.gz', '')
    )
