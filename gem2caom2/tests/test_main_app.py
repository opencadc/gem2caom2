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
import json
import logging
import pytest
import tempfile


from astropy.io.votable import parse_single_table

import gem2caom2.external_metadata as em

from gem2caom2 import main_app2, APPLICATION, ARCHIVE, SCHEME, svofps
from caom2.diff import get_differences
from caom2pipe import manage_composable as mc

from hashlib import md5
import os
import sys

from mock import patch

pytest.main(args=['-s', os.path.abspath(__file__)])
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')
INSTRUMENTS = ('GMOS', 'NIRI')
PLUGIN = os.path.join(os.path.dirname(THIS_DIR), 'main_app.py')

# structured by file id, observation id, filter_name (when looking up
# from SVOFPS for imaging files - x means there's no filter name), and
# instrument name
LOOKUP = {
    # GMOS
    'N20030107S0163': ['GN-2003A-Q-22-3-004', 'GMOS-N.i', 'GMOS'],
    'N20071219S0193': ['GN-2007B-Q-112-14-018', 'x', 'GMOS'],
    'N20090313S0180': ['GN-2009A-Q-21-115-001', 'GMOS-N.r', 'GMOS'],
    'N20100104S0208': ['GN-2009B-Q-121-15-001', 'x', 'GMOS'],
    'N20100104S0210': ['GN-2009B-Q-121-15-003', 'x', 'GMOS'],
    'N20100115S0346': ['GN-2010A-Q-35-10-002', 'x', 'GMOS'],
    'N20120105S0344': ['GN-2011A-Q-31-21-005', 'GMOS-N.g', 'GMOS'],
    'N20131203S0006': ['GN-2013B-Q-28-150-002', 'GMOS-N.g', 'GMOS'],
    'N20150217S0380': ['GN-2015A-C-2-96-002', 'GMOS-N.r', 'GMOS'],
    'N20150220S0320': ['GN-2015A-C-4-24-086', 'GMOS-N.r', 'GMOS'],
    'N20150216S0129': ['GN-2015A-Q-36-15-001', 'GMOS-N.i', 'GMOS'],
    'N20150216S0142': ['GN-2015A-Q-91-5-002', 'GMOS-N.r', 'GMOS'],
    'N20030104S0065': ['GN-CAL20030104-14-001', 'x', 'GMOS'],
    'N20030104S0161': ['GN-CAL20030104-18-003', 'GMOS-N.g', 'GMOS'],
    'N20150217S0274': ['GN-CAL20150217-2-003', 'x', 'GMOS'],
    'N20150929S0013': ['GN-CAL20150925-2-007', 'x', 'GMOS'],
    'S20040124S0077': ['GS-2003B-Q-5-29-003', 'x', 'GMOS'],
    'S20060101S0075': ['GS-2005B-Q-22-17-006', 'x', 'GMOS'],
    'S20060103S0143': ['GS-2005B-Q-22-29-004', 'x', 'GMOS'],
    'S20060128S0316': ['GS-2005B-Q-27-33-001', 'x', 'GMOS'],
    'S20060122S0004': ['GS-2005B-Q-27-4-001', 'x', 'GMOS'],
    'S20060131S0110': ['GS-2005B-Q-54-2-004', 'x', 'GMOS'],
    'S20060620S0066': ['GS-2006A-Q-48-15-002', 'x', 'GMOS'],
    'S20080526S0024': ['GS-2008A-Q-41-26-013', 'x', 'GMOS'],
    'S20090620S0145': ['GS-2009A-Q-30-6-007', 'x', 'GMOS'],
    'S20060125S0027': ['GS-CAL20060125-1-002', 'x', 'GMOS'],
    # NIRI
    'N20020620S0021': ['GN-2002A-C-5-1-001', 'x', 'NIRI'],
    'N20020620S0035': ['GN-2002A-C-5-1-015', 'x', 'NIRI'],
    'N20020620S0315': ['GN-2002A-C-5-21-002', 'x', 'NIRI'],
    'N20150404S0726': ['GN-2015A-C-1-20-001', 'x', 'NIRI'],
    'N20150404S0872': ['GN-2015A-C-1-27-001', 'x', 'NIRI'],
    'N20150405S0028': ['GN-2015A-C-1-27-071', 'x', 'NIRI'],
    # GSAIO
    'S20181023S0087': ['GS-CAL20181023-5-001', 'x', 'GSAIO'],
    # GNIRS
    'N20100915S0167': ['GN-2010B-Q-2-44-003', 'x', 'GNIRS'],
    'N20100722S0185': ['GN-2010B-SV-142-10-007', 'x', 'GNIRS'],
    'N20110323S0235': ['GN-2011A-Q-53-42-007', 'x', 'GNIRS'],
    'N20120104S0167': ['GN-2011B-Q-63-101-006', 'x', 'GNIRS'],
    'N20120102S0213': ['GN-2011B-Q-63-126-005', 'x', 'GNIRS'],
    'N20120101S0195': ['GN-2011B-Q-68-116-013', 'x', 'GNIRS'],
    'N20120103S0100': ['GN-2011B-Q-7-193-010', 'x', 'GNIRS'],
    'N20120103S0134': ['GN-2011B-Q-7-193-044', 'x', 'GNIRS'],
    'N20160123S0097': ['GN-2015B-SV-101-1061-005', 'x', 'GNIRS'],
    'N20151213S0022': ['GN-CAL20151213-6-002', 'x', 'GNIRS'],
    'N20160202S0098': ['GN-CAL20160202-3-039', 'x', 'GNIRS'],
    # GRACES
    'N20150807G0044': ['GN-2015B-Q-1-12-1003', 'x', 'GRACES'],
    'N20150604G0003': ['GN-CAL20150604-1000-1072', 'x', 'GRACES'],
    'N20150604G0014': ['GN-CAL20150604-1000-1081', 'x', 'GRACES'],
    'N20150807G0046': ['GN-CAL20150807-1000-1035', 'x', 'GRACES'],
}


def pytest_generate_tests(metafunc):
    if os.path.exists(TESTDATA_DIR):

        file_list = []
        # for root, dirs, files in os.walk(TESTDATA_DIR):
        for ii in ['GMOS', 'GNIRS', 'GRACES']:
        # for ii in ['GNIRS']:
            for root, dirs, files in os.walk('{}/{}'.format(TESTDATA_DIR, ii)):
                for file in files:
                    if file.endswith(".header"):
                        file_list.append(os.path.join(root, file))

        # metafunc.parametrize('test_name',
        # ['{}/GRACES/N20150807G0044.fits.header'.format(TESTDATA_DIR)])
        # metafunc.parametrize('test_name', file_list[8:])
        metafunc.parametrize('test_name', file_list)


def test_main_app(test_name):
    basename = os.path.basename(test_name)
    dirname = os.path.dirname(test_name)
    # file_id = basename.split('.fits')[0]
    file_id = _get_file_id(basename)
    product_id = LOOKUP[file_id][0]
    # lineage = mc.get_lineage(
    #     COLLECTION, product_id, '{}.fits'.format(file_id))
    lineage = _get_lineage(dirname, basename, product_id, file_id)
    input_file = '{}.in.xml'.format(product_id)
    actual_fqn = '{}/{}.actual.xml'.format(dirname, product_id)
    local = _get_local(test_name)
    plugin = PLUGIN

    with patch('caom2utils.fits2caom2.CadcDataClient') as data_client_mock, \
        patch('gem2caom2.external_metadata.get_obs_metadata') as gemini_client_mock, \
        patch('gem2caom2.svofps.get_votable') as svofps_mock:

        def get_file_info(archive, file_id):
            if '_prev' in file_id:
                return {'size': 10290,
                        'md5sum': 'md5:{}'.format(
                            md5('-37'.encode()).hexdigest()),
                        'type': 'image/jpeg'}
            else:
                return {'size': 665345,
                        'md5sum': 'md5:a347f2754ff2fd4b6209e7566637efad',
                        'type': 'application/fits'}

        def get_obs_metadata(obs_id):
            logging.error('here?')
            try:
                fname = '{}/{}/json/{}.json'.format(TESTDATA_DIR,
                                                    LOOKUP[file_id][2], obs_id)
                with open(fname) as f:
                    y = json.loads(f.read())
                    em.obs_metadata = y[0]
            except Exception as e:
                logging.error(e)

        def get_votable(url):
            x = url.split('/')
            filter_name = x[len(x) - 1]
            votable = parse_single_table(
                '{}/votable/{}.xml'.format(TESTDATA_DIR, filter_name))
            return votable, None

        logging.error('file_id is {}'.format(file_id))
        data_client_mock.return_value.get_file_info.side_effect = \
            get_file_info
        gemini_client_mock.return_value = get_obs_metadata(product_id)
        svofps_mock.return_value = get_votable(
            '{}/{}'.format(svofps.SVO_URL, LOOKUP[file_id][1]))

        sys.argv = \
            ('{} --no_validate --local {} '
             '--plugin {} --module {} --in {}/{} --out {} --lineage {}'.
             format(APPLICATION, local, plugin, plugin, dirname,
                    input_file, actual_fqn, lineage)).split()
        print(sys.argv)
        main_app2()
        expected_fqn = '{}/{}.xml'.format(dirname, product_id)
        expected = mc.read_obs_from_file(expected_fqn)
        actual = mc.read_obs_from_file(actual_fqn)
        result = get_differences(expected, actual, 'Observation')
        if result:
            msg = 'Differences found in observation {}\n{}'. \
                format(expected.observation_id, '\n'.join(
                [r for r in result]))
            raise AssertionError(msg)
        # assert False  # cause I want to see logging messages


def _build_temp_content(test_name):
    x = test_name.split('/')
    length = len(x)
    stuff = x[length-1].split('.')[0]
    temp_named_file = '/tmp/{}.fits.header'.format(stuff)
    content = None
    with open(test_name) as f:
        content = f.readlines()

    if content is not None:
        import caom2utils.fits2caom2 as f2c2
        with open(temp_named_file, 'w') as f:
            f.writelines(f2c2._make_understandable_string(content[0]))

    return temp_named_file


def _get_local(test_name):
    jpg = test_name.replace('.fits.header', '.jpg')
    # TODO - fix header metadata as returned by Gemini fullmetadata service
    # header_name = _build_temp_content(test_name)
    header_name = test_name
    if os.path.exists(jpg):
        return '{} {}'.format(jpg, header_name)
    else:
        return header_name


def _get_file_id(basename):
    if basename.endswith('jpg'):
        return basename.split('.jpg')[0]
    else:
        return basename.split('.fits')[0]


def _get_lineage(dirname, basename, product_id, file_id):
    jpg_file = basename.replace('.fits.header', '.jpg')
    if os.path.exists(os.path.join(dirname, jpg_file)):
        jpg = mc.get_lineage(ARCHIVE, product_id, '{}.jpg'.format(file_id),
                             SCHEME)
        fits = mc.get_lineage(ARCHIVE, product_id, '{}.fits'.format(file_id),
                              SCHEME)
        return '{} {}'.format(jpg, fits)
    else:
        return mc.get_lineage(ARCHIVE, product_id, '{}.fits'.format(file_id),
                              SCHEME)
