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

import os
import sys

import pytest

from shutil import copyfile

import gem2caom2.external_metadata as em

from gem2caom2 import main_app, builder
from gem2caom2.util import Inst
from caom2pipe import manage_composable as mc

from unittest.mock import patch, Mock
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
        ]:
            walk_dir = _get_inst_name(ii)
            for root, dirs, files in os.walk(
                f'{gem_mocks.TEST_DATA_DIR}/{walk_dir}'
            ):
                for file in files:
                    if file.endswith(".header"):
                        file_list.append(os.path.join(root, file))

        metafunc.parametrize('test_name', file_list)


@patch('gem2caom2.builder.defining_metadata_finder')
@patch('caom2utils.cadc_client_wrapper.StorageClientWrapper')
@patch('caom2utils.fits2caom2.Client')
@patch('caom2pipe.astro_composable.get_vo_table_session')
@patch('gem2caom2.program_metadata.get_pi_metadata')
@patch('gem2caom2.external_metadata.get_obs_metadata')
@patch('caom2pipe.client_composable.query_tap_client')
@patch('gem2caom2.external_metadata.CadcTapClient')
def test_main_app(
    client_mock,
    tap_mock,
    gemini_client_mock,
    gemini_pi_mock,
    svofps_mock,
    cadc_client_mock,
    get_file_info_mock,
    dmf_mock,
    test_name,
):
    # client_mock present because of global in external_metadata
    cadc_client_mock.get_node.side_effect = gem_mocks.mock_get_node
    gemini_client_mock.side_effect = gem_mocks.mock_get_obs_metadata
    gemini_pi_mock.side_effect = gem_mocks.mock_get_pi_metadata
    svofps_mock.side_effect = gem_mocks.mock_get_votable
    tap_mock.side_effect = gem_mocks.mock_query_tap
    get_file_info_mock.return_value.info.side_effect = (
        gem_mocks.mock_get_file_info
    )
    dmf_mock.get.side_effect = gem_mocks.mock_get_dm

    getcwd_orig = os.getcwd
    os.getcwd = Mock(
        return_value=os.path.join(gem_mocks.TEST_DATA_DIR, 'si_config'),
    )

    try:
        test_config = mc.Config()
        test_config.get_executors()
        test_config.features.supports_latest_client = True

        em.set_ofr(None)
        em.init_global(test_config)
        test_data_size = os.stat(
            os.path.join(gem_mocks.TEST_DATA_DIR, 'from_paul.txt')
        )
        app_size = os.stat('/app/data/from_paul.txt')
        if test_data_size.st_size != app_size.st_size:
            copyfile(
                os.path.join(gem_mocks.TEST_DATA_DIR, 'from_paul.txt'),
                '/app/data/from_paul.txt',
            )
        basename = os.path.basename(test_name)
        dirname = os.path.dirname(test_name)
        file_id = _get_file_id(basename)
        obs_id = _get_obs_id(file_id)
        product_id = file_id
        lineage = _get_lineage(dirname, basename, test_config)
        input_file = f'{product_id}.in.xml'
        actual_fqn = _get_actual_file_name(dirname, product_id)
        local = _get_local(test_name)
        plugin = gem_mocks.PLUGIN

        if os.path.exists(actual_fqn):
            os.remove(actual_fqn)

        if os.path.exists(os.path.join(dirname, input_file)):
            sys.argv = (
                f'{main_app.APPLICATION} --quiet --no_validate --local '
                f'{local} --plugin {plugin} --module {plugin} '
                f'--in {dirname}/{input_file} --out {actual_fqn} '
                f'--lineage {lineage} '
                f'--resource-id '
                f'{test_config.storage_inventory_resource_id}'
            ).split()
        else:
            sys.argv = (
                f'{main_app.APPLICATION} --quiet --no_validate --local '
                f'{local} --plugin {plugin} --module {plugin} '
                f'--observation {main_app.COLLECTION} {obs_id} '
                f'--out {actual_fqn} --lineage {lineage} '
                f'--resource-id '
                f'{test_config.storage_inventory_resource_id}'
            ).split()
        print(sys.argv)
        main_app.to_caom2()
        expected_fqn = _get_expected_file_name(dirname, product_id)

        compare_result = mc.compare_observations(actual_fqn, expected_fqn)

        if compare_result is not None:
            raise AssertionError(compare_result)
        # assert False  # cause I want to see logging messages
    finally:
        os.getcwd = getcwd_orig


def _get_obs_id(file_id):
    return gem_mocks.LOOKUP[file_id][0]


def _get_instr(file_id):
    temp = gem_mocks.LOOKUP[file_id][1]
    return _get_inst_name(temp)


def _get_program_id(file_id):
    return gem_mocks.LOOKUP[file_id][2]


def _get_local(test_name):
    jpg = test_name.replace('.fits.header', '.jpg')
    header_name = test_name
    if os.path.exists(jpg):
        return f'{jpg} {header_name}'
    else:
        return header_name


def _get_file_id(basename):
    if basename.endswith('jpg'):
        return basename.split('.jpg')[0]
    else:
        return basename.split('.fits')[0]


def _get_lineage(dirname, basename, config):
    jpg_file = basename.replace('.fits.header', '.jpg')
    name_builder = builder.GemObsIDBuilder(config)
    storage_name = name_builder.build(basename.replace('.header', ''))
    if os.path.exists(os.path.join(dirname, jpg_file)):
        jpg_storage_name = name_builder.build(jpg_file)
        jpg = jpg_storage_name.lineage
        fits = storage_name.lineage.replace('.header', '')
        result = f'{jpg} {fits}'
    else:
        result = storage_name.lineage.replace('.header', '')
    return result


def _get_expected_file_name(dirname, product_id):
    return f'{dirname}/{product_id}.expected.xml'


def _get_actual_file_name(dirname, product_id):
    return f'{dirname}/{product_id}.actual.xml'


def _get_inst_name(inst):
    walk_dir = inst
    if inst != 'processed' and isinstance(inst, Inst):
        walk_dir = inst.value
    return walk_dir
