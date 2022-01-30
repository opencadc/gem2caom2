# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
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
#  : 4 $
#
# ***********************************************************************
#

import logging

from tempfile import TemporaryDirectory
from cadcdata import FileInfo
from caom2.diff import get_differences
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from caom2pipe import reader_composable as rdc
from gem2caom2 import GemName, external_metadata, fits2caom2_augmentation

from mock import patch, Mock
import gem_mocks


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
):
    access_url.return_value = 'https://localhost:8080'
    get_pi_mock.side_effect = gem_mocks.mock_get_pi_metadata
    dmf_mock.return_value.get.side_effect = gem_mocks.mock_get_dm
    svofps_mock.side_effect = gem_mocks.mock_get_votable

    with TemporaryDirectory() as tmp_dir_name:

        test_f_name = 'gS20171114S0185_bias.fits.header'
        test_obs_id = 'GS-CAL20171114-2-086-G-BIAS'
        storage_name = GemName(obs_id=test_obs_id, file_name=test_f_name)
        test_fqn = f'{gem_mocks.TEST_DATA_DIR}/broken_files/{test_f_name}'
        storage_name.source_names = [test_fqn]

        test_config = mc.Config()
        test_config.task_types = [mc.TaskType.SCRAPE]
        test_config.use_local_files = True
        test_config.data_sources = [f'{gem_mocks.TEST_DATA_DIR}/broken_files']
        test_config.working_directory = tmp_dir_name
        test_config.proxy_fqn = f'{tmp_dir_name}/test_proxy.pem'

        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')

        external_metadata.get_gofr(test_config)
        observation = None
        expected_fqn = (
            f'{gem_mocks.TEST_DATA_DIR}/broken_files/{test_obs_id}.expected.xml'
        )
        file_info = FileInfo(
            id=storage_name.file_uri, file_type='application/fits'
        )
        headers = ac.make_headers_from_file(test_fqn)
        metadata_reader = rdc.FileMetadataReader()
        metadata_reader._headers = {storage_name.file_uri: headers}
        metadata_reader._file_info = {storage_name.file_uri: file_info}
        kwargs = {
            'storage_name': storage_name,
            'metadata_reader': metadata_reader,
        }
        logging.getLogger('caom2utils.fits2caom2').setLevel(logging.INFO)
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
