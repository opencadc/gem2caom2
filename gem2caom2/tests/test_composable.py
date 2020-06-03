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
import os
import pytest
import shutil
import sys

from datetime import datetime
from shutil import copyfile
from mock import patch, Mock
import gem_mocks

from caom2 import SimpleObservation, Algorithm
from caom2pipe import execute_composable as ec
from caom2pipe import manage_composable as mc
from gem2caom2 import composable, main_app, gem_name

STATE_FILE = '/usr/src/app/state.yml'
TODO_FILE = '/usr/src/app/todo.txt'
REJECTED_FILE = '/usr/src/app/logs/rejected.yml'
PROGRESS_FILE = '/usr/src/app/logs/progress.txt'
PUBLIC_TEST_JSON = f'{gem_mocks.TEST_DATA_DIR}/json/GN-2019B-ENG-1-160-008.json'


class MyExitError(Exception):
    pass


@pytest.fixture(scope='session', autouse=True)
def write_gemini_data_file():
    copyfile(os.path.join(gem_mocks.TEST_DATA_DIR, 'from_paul.txt'),
             '/app/data/from_paul.txt')


@patch('sys.exit', Mock(return_value=MyExitError))
def test_run_by_tap_query():
    # preconditions
    _write_state()

    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)

    try:
        # execution
        with patch('caom2pipe.manage_composable.query_tap') as \
                query_endpoint_mock, \
                patch('caom2pipe.execute_composable.run_by_file') \
                as run_mock:
            query_endpoint_mock.return_value = {
                'uri': ['gemini:GEM/2004may19_0255.fits']}
            composable.run_by_tap_query()
            assert run_mock.called, 'should have been called'
    finally:
        os.getcwd = getcwd_orig


@patch('caom2pipe.execute_composable.OrganizeExecutesWithDoOne.do_one')
def test_run(run_mock):
    test_obs_id = 'GS-2004A-Q-6-27-0255'
    test_f_id = '2004may19_0255'
    test_f_name = f'{test_f_id}.fits'
    _write_todo(test_f_name)

    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        # execution
        composable._run()
        assert run_mock.called, 'should have been called'
        args, kwargs = run_mock.call_args
        test_storage = args[0]
        assert isinstance(
            test_storage, gem_name.GemName), type(test_storage)
        assert test_storage.obs_id == test_obs_id, 'wrong obs id'
        assert test_storage.file_name == test_f_name, 'wrong file name'
        assert test_storage.fname_on_disk == test_f_name, \
            'wrong fname on disk'
        assert test_storage.url is None, 'wrong url'
        assert test_storage.lineage == \
            f'{test_f_id}/gemini:GEM/{test_f_name}', 'wrong lineage'
        assert test_storage.external_urls == \
            f'https://archive.gemini.edu/fullheader/{test_f_name}', \
            'wrong external urls'
    finally:
        os.getcwd = getcwd_orig


@patch('caom2pipe.execute_composable.OrganizeExecutesWithDoOne.do_one')
def test_run_errors(run_mock):
    test_obs_id = 'TX20131117_flt.3002'
    test_f_id = 'TX20131117_flt.3002'
    test_f_name = f'{test_f_id}.fits'
    _write_todo(test_f_name)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        composable._run()
        assert run_mock.called, 'should have been called'
        args, kwargs = run_mock.call_args
        test_storage = args[0]
        assert isinstance(
            test_storage, gem_name.GemName), type(test_storage)
        assert test_storage.obs_id == test_obs_id, 'wrong obs id'
        assert test_storage.file_name == test_f_name, 'wrong file name'
        assert test_storage.fname_on_disk == test_f_name, \
            'wrong fname on disk'
        assert test_storage.url is None, 'wrong url'
        assert test_storage.lineage == \
            f'{test_f_id}/gemini:GEM/{test_f_id}.fits', 'wrong lineage'
        assert test_storage.external_urls == \
            f'https://archive.gemini.edu/fullheader/{test_f_id}.fits', \
            'wrong external urls'
    finally:
        os.getcwd = getcwd_orig


@patch('caom2pipe.execute_composable.OrganizeExecutesWithDoOne.do_one')
@patch('caom2pipe.manage_composable.query_endpoint')
@patch('gem2caom2.external_metadata.get_obs_metadata')
@patch('caom2pipe.manage_composable.query_tap_client')
def test_run_incremental_rc(tap_mock, get_obs_mock, query_mock, run_mock):

    get_obs_mock.side_effect = gem_mocks.mock_get_obs_metadata
    query_mock.side_effect = gem_mocks.mock_query_endpoint_3
    tap_mock.side_effect = gem_mocks.mock_query_tap

    _write_state(prior_timestamp='2020-03-06 03:22:10.787835')
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        composable._run_by_incremental()
        assert run_mock.called, 'should have been called'
        args, kwargs = run_mock.call_args
        test_storage = args[0]
        assert isinstance(
            test_storage, gem_name.GemName), type(test_storage)
        import logging
        logging.error(test_storage)
        assert test_storage.obs_id == 'test_data_label', 'wrong obs id'
        assert test_storage.file_name == 'S20200303S0353.fits', \
            'wrong file_name'
        assert test_storage.file_id == 'S20200303S0353', 'wrong file_id'
        assert test_storage.fname_on_disk == 'S20200303S0353.fits', \
            'wrong fname on disk'
        assert test_storage.url is None, 'wrong url'
        # there are six files returned by the mock, and they each have the
        # same data label, so they all end up in this lineage
        assert test_storage.lineage == \
            'S20200303S0025/gemini:GEM/S20200303S0025.fits ' \
            'S20200303S0026/gemini:GEM/S20200303S0026.fits ' \
            'S20200303S0027/gemini:GEM/S20200303S0027.fits ' \
            'S20200303S0351/gemini:GEM/S20200303S0351.fits ' \
            'S20200303S0352/gemini:GEM/S20200303S0352.fits ' \
            'S20200303S0353/gemini:GEM/S20200303S0353.fits', 'wrong lineage'
        assert test_storage.external_urls == \
            'https://archive.gemini.edu/fullheader/S20200303S0025.fits ' \
            'https://archive.gemini.edu/fullheader/S20200303S0026.fits ' \
            'https://archive.gemini.edu/fullheader/S20200303S0027.fits ' \
            'https://archive.gemini.edu/fullheader/S20200303S0351.fits ' \
            'https://archive.gemini.edu/fullheader/S20200303S0352.fits ' \
            'https://archive.gemini.edu/fullheader/S20200303S0353.fits', \
            'wrong external urls'
    except Exception as e:
        assert False, e
    finally:
        os.getcwd = getcwd_orig


@patch('sys.exit', Mock(return_value=MyExitError))
def test_run_by_tap_query_2():
    test_obs_id = 'GS-2017A-Q-58-66-027'
    test_f_id = 'S20170905S0318'
    test_f_name = f'{test_f_id}.fits'
    _write_state()
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        # execution
        with patch('caom2pipe.execute_composable._do_one') as run_mock, \
                patch('gem2caom2.work.TapNoPreviewQuery') as query_mock:
            query_mock.return_value.todo.return_value = [test_f_name]
            composable.run_by_tap_query()
            assert run_mock.called, 'should have been called'
            args, kwargs = run_mock.call_args
            assert args[3] == 'gem2caom2', 'wrong command'
            test_storage = args[2]
            assert isinstance(test_storage, gem_name.GemName), \
                type(test_storage)
            assert test_storage.obs_id == test_obs_id, 'wrong obs id'
            assert test_storage.file_name == test_f_name, 'wrong file name'
            assert test_storage.fname_on_disk == test_f_name, \
                'wrong fname on disk'
            assert test_storage.url is None, 'wrong url'
            assert test_storage.lineage == \
                f'{test_f_id}/gemini:GEM/{test_f_name}', 'wrong lineage'
            assert test_storage.external_urls == \
                'https://archive.gemini.edu/fullheader/{}.fits'.format(
                       test_f_id), 'wrong external urls'

            args, kwargs = query_mock.call_args
            test_config = args[1]
            assert isinstance(test_config, mc.Config), 'wrong arg type'

            args, kwargs = query_mock.return_value.todo.call_args
            test_arg = args[0]
            assert isinstance(test_arg, datetime), type(test_arg)
    finally:
        os.getcwd = getcwd_orig


@patch('sys.exit', Mock(return_value=MyExitError))
def test_run_by_tap_query_rejected_bad_metadata():
    test_obs_id = 'GS-2017A-Q-58-66-027'
    _write_state()
    _write_rejected(test_obs_id)

    if os.path.exists(PROGRESS_FILE):
        os.unlink(PROGRESS_FILE)

    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        # execution
        with patch('gem2caom2.work.TapNoPreviewQuery') as query_mock:
            query_mock.return_value.todo.return_value = [test_obs_id]
            composable.run_by_tap_query()
            args, kwargs = query_mock.return_value.todo.call_args
            test_time = args[0]
            assert isinstance(test_time, datetime), type(test_time)
            assert os.path.exists(PROGRESS_FILE), 'should log'
            args, kwargs = query_mock.call_args
            test_config = args[1]
            assert isinstance(test_config, mc.Config), type(test_config)
            assert test_config.state_fqn == STATE_FILE, 'wrong state file'
    finally:
        os.getcwd = getcwd_orig
        if os.path.exists(REJECTED_FILE):
            os.unlink(REJECTED_FILE)


def test_run_by_in_memory_query():
    _write_state('2018-12-19T20:53:16')
    test_obs_id = 'GN-2015A-Q-36-15-001'
    test_f_id = 'N20150216S0129'
    if os.path.exists(TODO_FILE):
        os.unlink(TODO_FILE)

    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        # execution
        with patch('caom2pipe.execute_composable._do_one') \
                as run_mock:
            sys.argv = ['mock']
            composable._run_by_in_memory()
            assert run_mock.called, 'should have been called'
            args, kwargs = run_mock.call_args
            assert args[3] == main_app.APPLICATION, 'wrong command'
            test_storage = args[2]
            assert isinstance(test_storage, gem_name.GemName), \
                type(test_storage)
            assert test_storage.obs_id == test_obs_id, 'wrong obs id'
            assert test_storage.file_name == 'N20150216S0129.fits', \
                'wrong file name'
            assert test_storage.fname_on_disk == 'N20150216S0129.fits', \
                'wrong fname on disk'
            assert test_storage.url is None, 'wrong url'
            assert test_storage.lineage == \
                   '{}/gemini:GEM/{}.fits'.format(test_f_id, test_f_id), \
                'wrong lineage'
            assert test_storage.external_urls == \
                   'https://archive.gemini.edu/fullheader/{}.fits'.format(
                       test_f_id), 'wrong external urls'
    finally:
        os.getcwd = getcwd_orig


@patch('sys.exit', Mock(return_value=MyExitError))
@patch('caom2pipe.manage_composable.exec_cmd')
@patch('caom2pipe.execute_composable.CAOM2RepoClient')
@patch('caom2pipe.execute_composable.CadcDataClient')
@patch('caom2pipe.manage_composable.read_obs_from_file')
@patch('caom2pipe.manage_composable.query_endpoint')
@pytest.mark.skip('waiting for gemini incremental endpoint')
def test_run_by_incremental2(query_mock, read_mock,
                             data_client_mock, repo_mock, exec_mock):
    data_client_mock.return_value.get_file_info.side_effect = \
        gem_mocks.mock_get_file_info
    data_client_mock.return_value.get_file.side_effect = Mock()
    exec_mock.side_effect = Mock()
    repo_mock.return_value.create.side_effect = gem_mocks.mock_repo_create
    repo_mock.return_value.read.side_effect = gem_mocks.mock_repo_read
    repo_mock.return_value.update.side_effect = gem_mocks.mock_repo_update

    def _read_mock(ignore_fqn):
        return SimpleObservation(collection='TEST',
                                 observation_id='TEST_OBS_ID',
                                 algorithm=Algorithm('exposure'))
    read_mock.side_effect = _read_mock

    def _query_mock():
        x = json.loads("""[{
        "name": "2002feb11_0180.fits",
        "filename": "2002feb11_0180.fits.bz2",
        "path": "",
        "compressed": true,
        "file_size": 397169,
        "data_size": 1056960,
        "file_md5": "850d2df8b71c894984b4ff9d89cebcd8",
        "data_md5": "cf3ccf17539b23c4aa68c80308f885e2",
        "lastmod": "2019-12-17 00:21:08.159934+00:00",
        "mdready": true,
        "entrytime": "2019-12-17 00:21:08.178516+00:00",
        "size": 397169,
        "md5": "850d2df8b71c894984b4ff9d89cebcd8",
        "program_id": "GS-2002A-Q-8",
        "engineering": false,
        "science_verification": false,
        "calibration_program": false,
        "observation_id": "GS-2002A-Q-8-4",
        "data_label": "GS-2002A-Q-8-4-0180",
        "telescope": "Gemini-South",
        "instrument": "PHOENIX",
        "ut_datetime": "2002-04-04 21:47:53",
        "local_time": null,
        "observation_type": null,
        "observation_class": null,
        "object": "SZ Cha",
        "ra": 164.574204166667,
        "dec": -77.287963888889,
        "azimuth": null,
        "elevation": null,
        "cass_rotator_pa": null,
        "airmass": 1.478,
        "filter_name": "2150_(2)",
        "exposure_time": 240.0,
        "disperser": null,
        "camera": null,
        "central_wavelength": null,
        "wavelength_band": null,
        "focal_plane_mask": "107u_1.0-5.0 (8)",
        "detector_binning": "1x1",
        "detector_gain_setting": null,
        "detector_roi_setting": "Fixed",
        "detector_readspeed_setting": null,
        "detector_welldepth_setting": null,
        "detector_readmode_setting": "None",
        "spectroscopy": true,
        "mode": "spectroscopy",
        "adaptive_optics": false,
        "laser_guide_star": false,
        "wavefront_sensor": null,
        "gcal_lamp": null,
        "raw_iq": 50,
        "raw_cc": 50,
        "raw_wv": 80,
        "raw_bg": 80,
        "requested_iq": null,
        "requested_cc": null,
        "requested_wv": null,
        "requested_bg": null,
        "qa_state": "Undefined",
        "release": "2003-08-11",
        "reduction": "RAW",
        "types": "{'UNPREPARED', 'SOUTH', 'RAW', 'PHOENIX', 'GEMINI', 'SPECT'}",
        "phot_standard": null,
        "results_truncated": true
    }]""")
        return x

    query_result = gem_mocks.Object()
    query_result.json = _query_mock
    query_mock.return_value = query_result

    _write_cert()
    prior_s = datetime.utcnow().timestamp() - 60
    _write_state(prior_s)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=f'{gem_mocks.TEST_DATA_DIR}/edu_query')
    try:
        # execution
        test_result = composable._run_by_incremental()
        assert test_result == 0, 'wrong result'
    finally:
        os.getcwd = getcwd_orig

    assert repo_mock.return_value.create.called, 'create not called'
    assert repo_mock.return_value.read.called, 'read not called'
    assert exec_mock.called, 'exec mock not called'
    param, level_as = ec.CaomExecute._specify_logging_level_param(
        logging.ERROR)
    py_version = f'{sys.version_info.major}.{sys.version_info.minor}'
    exec_mock.assert_called_with(
        (f'gem2caom2 --quiet --cert /usr/src/app/cadcproxy.pem '
         f'--observation GEMINI GS-2002A-Q-8-4-0180 '
         f'--out /usr/src/app/logs/GS-2002A-Q-8-4-0180.fits.xml '
         f'--external_url '
         f'https://archive.gemini.edu/fullheader/2002feb11_0180.fits '
         f'--plugin '
         f'/usr/local/lib/python{py_version}/site-packages/gem2caom2/'
         f'gem2caom2.py '
         f'--module '
         f'/usr/local/lib/python{py_version}/site-packages/gem2caom2/'
         f'gem2caom2.py '
         f'--lineage 2002feb11_0180/gemini:GEM/2002feb11_0180.fits'),
        level_as), \
    'exec mock wrong parameters'
    assert not data_client_mock.return_value.get_file_info.called, \
        'data client mock get file info should not be not called'
    assert query_mock.called, 'query mock not called'


@patch('caom2pipe.execute_composable.OrganizeExecutesWithDoOne.do_one')
@patch('caom2pipe.manage_composable.query_tap_client')
def test_run_by_rc_public(tap_mock, exec_mock):
    exec_mock.side_effect = Mock(return_value=0)
    tap_mock.side_effect = gem_mocks.mock_query_tap
    expected_fqn = f'/usr/src/app/logs/{gem_mocks.TEST_BUILDER_OBS_ID}' \
                   f'.expected.xml'
    if not os.path.exists(expected_fqn):
        shutil.copy(
            f'{gem_mocks.TEST_DATA_DIR}/expected.xml', expected_fqn)

    _write_cert()
    prior_s = datetime.utcnow().timestamp() - 1440 * 60
    _write_state(prior_s)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=f'{gem_mocks.TEST_DATA_DIR}/edu_query')
    test_f_id = 'N20191101S0007'
    try:
        # execution
        sys.argv = ['test command']
        test_result = composable._run_rc_state_public()
        assert test_result == 0, 'wrong result'
    except Exception as e:
        import logging
        import traceback
        logging.error(traceback.format_exc())
    finally:
        os.getcwd = getcwd_orig

    assert exec_mock.called, 'exec mock not called'
    args, kwargs = exec_mock.call_args
    test_storage = args[0]
    assert isinstance(
        test_storage, gem_name.GemName), type(test_storage)
    assert test_storage.obs_id == 'GN-2019B-ENG-1-160-008', 'wrong obs id'
    assert test_storage.file_name == 'N20191101S0007.fits', 'wrong file_name'
    assert test_storage.file_id == test_f_id, 'wrong file_id'
    assert test_storage.fname_on_disk == f'{test_f_id}.fits', \
        'wrong fname on disk'
    assert test_storage.url is None, 'wrong url'
    assert test_storage.lineage == \
        f'{test_f_id}/gemini:GEM/{test_f_id}.fits', 'wrong lineage'
    assert test_storage.external_urls == \
        f'https://archive.gemini.edu/fullheader/{test_f_id}.fits', \
        'wrong external urls'
    assert tap_mock.called, 'tap mock not called'


def _write_todo(test_id):
    with open(TODO_FILE, 'w') as f:
        f.write('{}\n'.format(test_id))


def _write_state(prior_timestamp=None):
    # to ensure at least one spin through the execution loop, test case
    # must have a starting time greater than one config.interval prior
    # to 'now', default interval is 10 minutes
    if prior_timestamp is None:
        prior_s = datetime.utcnow().timestamp() - 15 * 60
    else:
        if type(prior_timestamp) is float:
            prior_s = prior_timestamp
        else:
            prior_s = mc.make_seconds(prior_timestamp)
    test_start_time = datetime.fromtimestamp(prior_s).isoformat()
    logging.error(f'test_start_time {test_start_time}')
    test_bookmark = {'bookmarks':
                         {'gemini_timestamp':
                              {'last_record': test_start_time}
                          }
                     }
    mc.write_as_yaml(test_bookmark, STATE_FILE)


def _write_rejected(test_obs_id):
    content = {'bad_metadata': [test_obs_id]}
    mc.write_as_yaml(content, REJECTED_FILE)


def _write_cert():
    if not os.path.exists('/usr/src/app/cadcproxy.pem'):
        with open('/usr/src/app/cadcproxy.pem', 'w') as f:
            f.write('cadc proxy content')
