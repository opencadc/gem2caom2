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

import os
import pytest
import shutil

from astropy.io import fits
from astropy.table import Table

from mock import patch, Mock

from caom2pipe import manage_composable as mc
from gem2caom2 import external_metadata as ext_md
from gem2caom2.util import Inst
from gem2caom2.obs_metadata import json_lookup


import gem_mocks


@pytest.mark.skip('')
@patch('gem2caom2.external_metadata.get_obs_metadata')
@patch('caom2pipe.client_composable.query_tap_client')
def test_caching_relationship(tap_mock, get_obs_mock):
    shutil.copyfile(
        f'{gem_mocks.TEST_DATA_DIR}/from_paul.txt', '/app/data/from_paul.txt'
    )
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        test_config = mc.Config()
        test_config.get_executors()
        ext_md.init_global(config=test_config)
        initial_length = 525
        tap_mock.side_effect = gem_mocks._query_mock_none
        get_obs_mock.side_effect = gem_mocks.mock_get_obs_metadata
        test_subject = ext_md.CachingObsFileRelationship()
        test_subject.tap_client = Mock()
        # test an entry that's not in the file, not at CADC, is at
        # archive.gemini.edu
        assert (
            len(test_subject.name_list) == initial_length
        ), 'bad initial length'
        test_result = test_subject.get_obs_id('N20200210S0077')
        assert test_result is not None, 'expect a gemini result'
        assert test_result == 'GN-CAL20200210-22-076', 'wrong gemini result'
        assert (
            len(test_subject.name_list) == initial_length + 1
        ), 'bad updated length from Gemini'

        # entry is not in file, but is at CADC
        tap_mock.side_effect = gem_mocks.mock_query_tap
        test_result = test_subject.get_obs_id('x')
        assert test_result is not None, 'expect a cadc result'
        assert test_result == 'test_data_label', 'wrong cadc result'
        assert (
            len(test_subject.name_list) == initial_length + 2
        ), 'bad updated length from cadc'

        # entry is in file
        test_result = test_subject.get_obs_id('N20170616S0540')
        assert test_result is not None, 'expect a file result'
        assert test_result == 'GN-CAL20170616-11-022', 'wrong file result'
        assert (
            len(test_subject.name_list) == initial_length + 2
        ), 'bad updated length from file'
    finally:
        os.getcwd = getcwd_orig


@pytest.mark.skip('')
def test_caching_relationship_unconnected():
    test_config = mc.Config()
    test_config.use_local_files = True
    test_config.task_types = [mc.TaskType.SCRAPE]
    test_subject = ext_md.CachingObsFileRelationship(test_config)
    # test an entry that's a file header on disk
    test_result = test_subject.get_obs_id('get_obs_id_from_file_on_disk')
    assert test_result is not None, 'expected result'
    assert test_result == 'GN-2006A-Q-90-1-001-MRG-ADD', 'wrong result'


@pytest.mark.skip('')
@patch('requests.Session')
@patch('gem2caom2.external_metadata.CadcTapClient')
def test_get_obs_metadata_not_at_gemini(tap_client_mock, session_mock):
    session_mock.get.side_effect = gem_mocks.mock_session_get_not_found
    test_config = mc.Config()
    test_config.working_directory = gem_mocks.TEST_DATA_DIR
    test_config.proxy_file_name = 'test_proxy.pem'
    ext_md.init_global(config=test_config)
    with pytest.raises(
        mc.CadcException, match=f'Could not find JSON record *'
    ):
        test_result = ext_md.get_obs_metadata('test_file_id')


@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('caom2utils.cadc_client_wrapper.get_local_file_headers', autospec=True)
@patch('caom2pipe.client_composable.query_tap_client', autospec=True)
@patch('gem2caom2.external_metadata.get_obs_metadata', autospec=True)
def test_dm_finder(get_obs_mock, caom2_mock, local_mock, cap_mock):
    cap_mock.return_value = 'https://localhost'
    test_file_id = 'rN20123456S9876'
    test_uri = f'gemini:GEMINI/{test_file_id}.fits'
    repaired_data_label = 'GN-2012A-B-012-345-6'
    test_data_label = f'{repaired_data_label}-R'
    json_lookup.flush()

    def _get_obs_md_mock(ignore):
        md = [
            {
                'data_label': test_data_label,
                'filename': f'{test_file_id}.fits',
                'lastmod': '2020-02-25T20:36:31.230',
                'instrument': 'GMOS',
            },
        ]
        json_lookup.add(md, test_file_id)
    get_obs_mock.side_effect = _get_obs_md_mock

    def _caom2_mock(ignore1, ignore2):
        return Table.read(
            f'observationID,instrument_name\n'
            f'{test_data_label},'
            f'GMOS\n'.split('\n'),
            format='csv',
        )
    caom2_mock.side_effect = _caom2_mock

    def _local_mock(ignore):
        hdr = fits.Header()
        hdr['DATALAB'] = test_data_label
        hdr['INSTRUME'] = 'GMOS'
        return [hdr]
    local_mock.side_effect = _local_mock
    os_path_exists_orig = os.path.exists
    os.path.exists = Mock(return_value=True)

    test_config = mc.Config()
    test_config.data_sources = [gem_mocks.TEST_DATA_DIR]
    test_config.proxy_fqn = os.path.join(
        gem_mocks.TEST_DATA_DIR, 'cadcproxy.pem'
    )
    test_config.tap_id = 'ivo://cadc.nrc.ca/test'

    try:
        for test_use_local in [True, False]:
            for test_connected in [True, False]:
                test_config.use_local_files = test_use_local
                if test_connected:
                    test_config.task_types = [mc.TaskType.VISIT]
                else:
                    test_config.task_types = [mc.TaskType.SCRAPE]

                test_subject = ext_md.DefiningMetadataFinder(test_config)
                assert test_subject is not None, (
                    f'ctor does not work:: '
                    f'local {test_use_local}, '
                    f'connected {test_connected}'
                )
                test_result = test_subject.get(test_uri)
                assert test_result is not None, (
                    f'expect a result '
                    f'local {test_use_local}, '
                    f'connected {test_connected}'
                )
                assert test_result.instrument is Inst.GMOS, (
                    f'instrument should be GMOS '
                    f'local {test_use_local}, '
                    f'connected {test_connected}'
                )
                assert test_result.data_label == repaired_data_label, (
                    f'data_label should be {repaired_data_label} '
                    f'local {test_use_local}, '
                    f'connected {test_connected}'
                )
    finally:
        os.path.exists = os_path_exists_orig
