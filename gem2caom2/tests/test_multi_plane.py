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
import json
import logging
import pytest


from astropy.io.votable import parse_single_table

import gem2caom2.external_metadata as em

from gem2caom2 import main_app2, APPLICATION, ARCHIVE, SCHEME
from caom2.diff import get_differences
from caom2pipe import manage_composable as mc


from hashlib import md5
import os
import sys

from mock import patch

pytest.main(args=['-s', os.path.abspath(__file__)])
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
PLUGIN = os.path.join(os.path.dirname(THIS_DIR), 'main_app.py')

# structured by observation id, list of file ids that make up a multi-plane
# observation
DIR_NAME = 'multi_plane'
LOOKUP = {'GS-CAL20101028-5-004': ['mrgS20101028S0134', 'S20101028S0134'],
          'GN-2007B-C-6-5-005': ['TX20071021_RAW.2037', 'TX20071021_SUM.2037'],
          'TX20170321.2505': ['TX20170321_raw.2505', 'TX20170321_red.2505',
                              'TX20170321_sum.2505'],
          'TX20170321.2507': ['TX20170321_raw.2507', 'TX20170321_red.2507'],
          'GS-2002B-Q-22-13-0161': ['2002dec02_0161',
                                    'P2002DEC02_0161_SUB.0001',
                                    'P2002DEC02_0161_SUB'],
          'GN-2015B-Q-1-12-1003': ['N20150807G0044', 'N20150807G0044i',
                                   'N20150807G0044m']
          }


def pytest_generate_tests(metafunc):
    obs_id_list = []
    for ii in LOOKUP:
        obs_id_list.append(ii)
    metafunc.parametrize('test_name', obs_id_list)


def test_multi_plane(test_name):
    obs_id = test_name
    lineage = _get_lineage(obs_id)
    input_file = '{}/{}/{}.in.xml'.format(TEST_DATA_DIR, DIR_NAME, obs_id)
    actual_fqn = '{}/{}/{}.actual.xml'.format(TEST_DATA_DIR, DIR_NAME, obs_id)

    local = _get_local(test_name)
    plugin = PLUGIN

    with patch('caom2utils.fits2caom2.CadcDataClient') as data_client_mock, \
            patch('gem2caom2.external_metadata.get_obs_metadata') as gemini_client_mock, \
            patch('gem2caom2.external_metadata.get_pi_metadata') as gemini_pi_mock, \
            patch('gem2caom2.svofps.get_vo_table') as svofps_mock:

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

        def get_obs_metadata(file_id):
            try:
                logging.error('obs metadata file_id {}'.format(file_id))
                fname = '{}/{}/json/{}.json'.format(TEST_DATA_DIR, DIR_NAME,
                                                    file_id)
                with open(fname) as f:
                    y = json.loads(f.read())
                    em.om.add(y, file_id)
            except Exception as e:
                logging.error(e)
                import traceback
                tb = traceback.format_exc()
                logging.error(tb)

        def get_pi_metadata(program_id):
            try:
                logging.error('program id is {}'.format(program_id))
                fname = '{}/{}/program/{}.xml'.format(TEST_DATA_DIR, DIR_NAME,
                                                      program_id)
                with open(fname) as f:
                    y = f.read()
                    from bs4 import BeautifulSoup
                    soup = BeautifulSoup(y, 'lxml')
                    tds = soup.find_all('td')
                    if len(tds) > 0:
                        title = tds[1].contents[0].replace('\n', ' ')
                        pi_name = tds[3].contents[0]
                        metadata = {'title': title,
                                    'pi_name': pi_name}
                        return metadata
                return None
            except Exception as e:
                logging.error(e)
                import traceback
                tb = traceback.format_exc()
                logging.error(tb)

        def mock_get_votable(url):
            try:
                x = url.split('/')
                filter_name = x[-1].replace('&VERB=0', '')
                votable = parse_single_table(
                    '{}/votable/{}.xml'.format(TEST_DATA_DIR, filter_name))
                return votable, None
            except Exception as e:
                logging.error('get_vo_table failure for url {}'.format(url))
                logging.error(e)
                return None, None

        data_client_mock.return_value.get_file_info.side_effect = \
            get_file_info
        gemini_client_mock.side_effect = get_obs_metadata
        gemini_pi_mock.side_effect = get_pi_metadata
        svofps_mock.side_effect = mock_get_votable

        if os.path.exists(actual_fqn):
            os.remove(actual_fqn)

        logging.error('{}'.format(os.path.join(DIR_NAME, input_file)))
        sys.argv = \
            ('{} --quiet --no_validate --local {} '
             '--plugin {} --module {} --in {} --out {} --lineage {}'.
             format(APPLICATION, local, plugin, plugin,
                    input_file, actual_fqn, lineage)).split()
        print(sys.argv)
        main_app2()
        expected_fqn = '{}/{}/{}.xml'.format(TEST_DATA_DIR, DIR_NAME, obs_id)
        expected = mc.read_obs_from_file(expected_fqn)
        actual = mc.read_obs_from_file(actual_fqn)
        result = get_differences(expected, actual, 'Observation')
        if result:
            msg = 'Differences found obs id {} n{}'.format(
                expected.observation_id, '\n'.join([r for r in result]))
            raise AssertionError(msg)
        # assert False  # cause I want to see logging messages


def _get_lineage(obs_id):
    result = ''
    for ii in LOOKUP[obs_id]:
        fits = mc.get_lineage(ARCHIVE, ii, '{}.fits'.format(ii), SCHEME)
        result = '{} {}'.format(result, fits)
    return result


def _get_local(obs_id):
    result = ''
    for ii in LOOKUP[obs_id]:
        result = '{} {}/{}/{}.fits.header'.format(result, TEST_DATA_DIR,
                                                  DIR_NAME, ii)
    return result
