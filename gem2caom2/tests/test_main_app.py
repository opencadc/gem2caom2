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

import logging
import os

from tempfile import TemporaryDirectory

import pytest

import gem2caom2.external_metadata as em

from cadcdata import FileInfo
from caom2.diff import get_differences
from gem2caom2 import builder, fits2caom2_augmentation, gemini_reader
from gem2caom2.util import Inst
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from caom2pipe import reader_composable as rdc

from unittest.mock import patch
import gem_mocks


pytest.main(args=['-s', os.path.abspath(__file__)])


def pytest_generate_tests(metafunc):
    if os.path.exists(gem_mocks.TEST_DATA_DIR):

        file_list = []
        for ii in [
            Inst.GMOS,
            Inst.NIRI,
            Inst.GPI,
            Inst.F2,
            Inst.GSAOI,
            Inst.NICI,
            Inst.TRECS,
            Inst.MICHELLE,
            Inst.GRACES,
            Inst.NIFS,
            Inst.GNIRS,
            Inst.PHOENIX,
            Inst.FLAMINGOS,
            Inst.HRWFS,
            Inst.HOKUPAA,
            Inst.OSCIR,
            Inst.BHROS,
            Inst.CIRPASS,
            Inst.TEXES,
            'processed',
            Inst.ALOPEKE,
            Inst.ZORRO,
            Inst.IGRINS,
        ]:
            walk_dir = _get_inst_name(ii)
            for root, dirs, files in os.walk(
                f'{gem_mocks.TEST_DATA_DIR}/{walk_dir}'
            ):
                for file in files:
                    if file.endswith(".header"):
                        file_list.append(os.path.join(root, file))

        metafunc.parametrize('test_name', file_list)


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

        em.get_gofr(test_config)
        test_reader = rdc.FileMetadataReader()
        test_metadata = gemini_reader.GeminiMetadataLookup(test_reader)
        test_builder = builder.GemObsIDBuilder(
            test_config, test_reader, test_metadata
        )
        storage_name = test_builder.build(test_name)
        file_info = FileInfo(
            id=storage_name.file_uri, file_type='application/fits'
        )
        headers = ac.make_headers_from_file(test_name)
        metadata_reader = rdc.FileMetadataReader()
        metadata_reader._headers = {storage_name.file_uri: headers}
        metadata_reader._file_info = {storage_name.file_uri: file_info}
        kwargs = {
            'storage_name': storage_name,
            'metadata_reader': metadata_reader,
        }
        logging.getLogger('caom2utils.fits2caom2').setLevel(logging.INFO)
        logging.getLogger('root').setLevel(logging.INFO)
        observation = None
        expected_fqn = (
            f'{os.path.dirname(test_name)}/{storage_name.file_id}.expected.xml'
        )
        in_fqn = expected_fqn.replace('.expected.', '.in.')
        if os.path.exists(in_fqn):
            observation = mc.read_obs_from_file(in_fqn)
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


# def _get_obs_id(file_id):
#     return gem_mocks.LOOKUP[file_id][0]


# def _get_instr(file_id):
#     temp = gem_mocks.LOOKUP[file_id][1]
#     return _get_inst_name(temp)


# def _get_program_id(file_id):
#     return gem_mocks.LOOKUP[file_id][2]


# def _get_local(test_name):
#     jpg = test_name.replace('.fits.header', '.jpg')
#     header_name = test_name
#     if os.path.exists(jpg):
#         return f'{jpg} {header_name}'
#     else:
#         return header_name


# def _get_file_id(basename):
#     if basename.endswith('jpg'):
#         return basename.split('.jpg')[0]
#     else:
#         return basename.split('.fits')[0]


# def _get_lineage(dirname, basename, config):
#     jpg_file = basename.replace('.fits.header', '.jpg')
#     name_builder = builder.GemObsIDBuilder(config)
#     storage_name = name_builder.build(basename.replace('.header', ''))
#     if os.path.exists(os.path.join(dirname, jpg_file)):
#         jpg_storage_name = name_builder.build(jpg_file)
#         jpg = jpg_storage_name.lineage
#         fits = storage_name.lineage.replace('.header', '')
#         result = f'{jpg} {fits}'
#     else:
#         result = storage_name.lineage.replace('.header', '')
#     return result


# def _get_expected_file_name(dirname, product_id):
#     return f'{dirname}/{product_id}.expected.xml'


# def _get_actual_file_name(dirname, product_id):
#     return f'{dirname}/{product_id}.actual.xml'
#
#
def _get_inst_name(inst):
    walk_dir = inst
    if inst != 'processed' and isinstance(inst, Inst):
        walk_dir = inst.value
    return walk_dir
