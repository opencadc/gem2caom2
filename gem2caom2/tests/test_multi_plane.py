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

import logging
import os
import pytest

from tempfile import TemporaryDirectory
from caom2.diff import get_differences
from cadcdata import FileInfo
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from caom2pipe import reader_composable as rdc
from gem2caom2 import external_metadata, builder, gemini_reader
from gem2caom2 import fits2caom2_augmentation

from unittest.mock import patch, Mock

import gem_mocks

pytest.main(args=['-s', os.path.abspath(__file__)])

# structured by observation id, list of file ids that make up a multi-plane
# observation
DIR_NAME = 'multi_plane'
LOOKUP = {
    'GS-CAL20101028-5-004': ['mrgS20101028S0134', 'S20101028S0134'],
    'GN-2007B-C-6-5-005': ['TX20071021_RAW.2037', 'TX20071021_SUM.2037'],
    'TX20170321.2505': [
        'TX20170321_raw.2505',
        'TX20170321_red.2505',
        'TX20170321_sum.2505',
    ],
    'TX20170321.2507': ['TX20170321_raw.2507', 'TX20170321_red.2507'],
    'GS-2002B-Q-22-13-0161': [
        '2002dec02_0161',
        'P2002DEC02_0161_SUB.0001',
        'P2002DEC02_0161_SUB',
    ],
    'GN-2015B-Q-1-12-1003': [
        'N20150807G0044',
        'N20150807G0044i',
        'N20150807G0044m',
    ],
    'GN-2012A-Q-124-1-003': ['N20120905S0122', 'N20120905S0122_arc'],
    'GS-2006A-Q-60-11-001': ['S20060412S0056', 'rS20060412S0056'],
    'N20191219A0004': ['N20191219A0004b', 'N20191219A0004r'],
    'N20191219A0001': ['N20191219A0001b', 'N20191219A0001r'],
    'S20200316Z0308': ['S20200316Z0308b', 'S20200316Z0308r'],
    'S20190716Z0925': ['S20190716Z0925b', 'S20190716Z0925r'],
    'N20190313A0002': ['N20190313A0002b', 'N20190313A0002r'],
    'N20200215A0004': ['N20200215A0004b', 'N20200215A0004r'],
    'S20190912Z0264': ['S20190912Z0264b', 'S20190912Z0264r'],
    'GS-2003B-Q-23-17-001': ['S20050916S0159', 'rS20050916S0159'],
    'GS-2013B-Q-75-187-001': ['S20130922S0130', 'S20130922S0130_arc'],
    'GS-2020B-Q-211-52-1119': ['SDCH_20201119_0052', 'SDCK_20201119_0052'],
    'GS-CAL20200202-23-0201': ['SDCH_20200201_0023', 'SDCK_20200201_0023'],
    'GS-CAL20201123-20-1122': ['SDCH_20201122_0020', 'SDCK_20201122_0020'],
}


def pytest_generate_tests(metafunc):
    obs_id_list = []
    for ii in LOOKUP:
        obs_id_list.append(ii)
    metafunc.parametrize('test_name', obs_id_list)


@patch('caom2pipe.reader_composable.FileMetadataReader._retrieve_headers')
@patch('caom2pipe.astro_composable.get_vo_table_session')
@patch('gem2caom2.external_metadata.DefiningMetadataFinder')
@patch('gem2caom2.program_metadata.get_pi_metadata')
@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('gemProc2caom2.builder.CadcTapClient')
@patch('gem2caom2.external_metadata.CadcTapClient')
def test_visitor(
    em_tap_client_mock,
    builder_tap_client_mock,
    access_url,
    get_pi_mock,
    dmf_mock,
    svofps_mock,
    retrieve_headers_mock,
    test_name,
):
    access_url.return_value = 'https://localhost:8080'
    get_pi_mock.side_effect = gem_mocks.mock_get_pi_metadata
    dmf_mock.return_value.get.side_effect = gem_mocks.mock_get_dm
    svofps_mock.side_effect = gem_mocks.mock_get_votable
    retrieve_headers_mock.side_effect = gem_mocks._mock_headers

    with TemporaryDirectory() as tmp_dir_name:
        test_config = mc.Config()
        test_config.task_types = [mc.TaskType.SCRAPE]
        test_config.use_local_files = True
        test_config.data_sources = [os.path.dirname(test_name)]
        test_config.working_directory = tmp_dir_name
        test_config.proxy_fqn = f'{tmp_dir_name}/test_proxy.pem'

        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')

        external_metadata.get_gofr(test_config)
        observation = None
        expected_fqn = (
            f'{gem_mocks.TEST_DATA_DIR}/multi_plane/'
            f'{test_name}.expected.xml'
        )
        in_fqn = expected_fqn.replace('.expected.', '.in.')
        if os.path.exists(in_fqn):
            observation = mc.read_obs_from_file(in_fqn)
        for f_name in LOOKUP[test_name]:
            source_name = (
                f'{gem_mocks.TEST_DATA_DIR}/multi_plane/{f_name}.fits.header'
            )
            metadata_reader = rdc.FileMetadataReader()
            test_metadata = gemini_reader.GeminiMetadataLookup(metadata_reader)
            test_builder = builder.GemObsIDBuilder(
                test_config, metadata_reader, test_metadata
            )
            storage_name = test_builder.build(source_name)
            storage_name.source_names = [source_name]
            file_info = FileInfo(
                id=storage_name.file_uri, file_type='application/fits'
            )
            headers = ac.make_headers_from_file(source_name)
            metadata_reader._headers = {storage_name.file_uri: headers}
            metadata_reader._file_info = {storage_name.file_uri: file_info}
            kwargs = {
                'storage_name': storage_name,
                'metadata_reader': metadata_reader,
            }
            logging.getLogger('caom2utils.fits2caom2').setLevel(logging.INFO)
            logging.getLogger('TEXES').setLevel(logging.INFO)
            logging.getLogger('root').setLevel(logging.INFO)
            observation = fits2caom2_augmentation.visit(observation, **kwargs)

        expected = mc.read_obs_from_file(expected_fqn)
        compare_result = get_differences(expected, observation)
        if compare_result is not None:
            actual_fqn = expected_fqn.replace('expected', 'actual')
            mc.write_obs_to_file(observation, actual_fqn)
            compare_text = '\n'.join([r for r in compare_result])
            msg = (
                f'Differences found in observation {expected.observation_id}\n'
                f'{compare_text}. Check {actual_fqn}'
            )
            raise AssertionError(msg)
