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

from astropy.table import Column, Table
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
        with patch('caom2pipe.astro_composable.query_tap') as \
                query_endpoint_mock, \
                patch('caom2pipe.execute_composable.run_by_file') \
                as run_mock:
            query_endpoint_mock.return_value = {
                'observationID': ['GS-2004A-Q-6-27-0255']}
            composable.run_by_tap_query()
            assert run_mock.called, 'should have been called'
    finally:
        os.getcwd = getcwd_orig


@patch('sys.exit', Mock(return_value=MyExitError))
def test_run():
    test_obs_id = 'GS-2004A-Q-6-27-0255'
    test_f_id = '2004may19_0255'
    _write_todo(test_obs_id)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        # execution
        with patch('caom2pipe.execute_composable._do_one') \
                as run_mock:
            composable.run()
            assert run_mock.called, 'should have been called'
            args, kwargs = run_mock.call_args
            assert args[3] == 'gem2caom2', 'wrong command'
            test_storage = args[2]
            assert isinstance(
                test_storage, gem_name.GemName), type(test_storage)
            assert test_storage.obs_id == test_obs_id, 'wrong obs id'
            assert test_storage.file_name is None, 'wrong file name'
            assert test_storage.fname_on_disk is None, 'wrong fname on disk'
            assert test_storage.url is None, 'wrong url'
            assert test_storage.lineage == \
                '{}/gemini:GEM/{}.fits'.format(test_f_id, test_f_id), \
                'wrong lineage'
            assert test_storage.external_urls == \
                   'https://archive.gemini.edu/fullheader/{}.fits'.format(
                       test_f_id), 'wrong external urls'
    finally:
        os.getcwd = getcwd_orig


# gem2caom2 --verbose --cert /usr/src/app/cadcproxy.pem --observation GEMINI
# TX20131117_flt.3002 --out /usr/src/app/logs/TX20131117_flt.3002.fits.xml
# --plugin /usr/local/lib/python3.6/site-packages/gem2caom2/gem2caom2.py
# --module /usr/local/lib/python3.6/site-packages/gem2caom2/gem2caom2.py
# --lineage TX20131117_flt.3002/gemini:GEM/TX20131117_flt.3002.fits
@patch('sys.exit', Mock(return_value=MyExitError))
def test_run_errors():
    test_obs_id = 'TX20131117_flt.3002'
    test_f_id = 'TX20131117_flt.3002'
    _write_todo(test_obs_id)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        # execution
        with patch('caom2pipe.execute_composable._do_one') \
                as run_mock:
            composable.run()
            assert run_mock.called, 'should have been called'
            args, kwargs = run_mock.call_args
            assert args[3] == 'gem2caom2', 'wrong command'
            test_storage = args[2]
            assert isinstance(
                test_storage, gem_name.GemName), type(test_storage)
            assert test_storage.obs_id == test_obs_id, 'wrong obs id'
            assert test_storage.file_name is None, 'wrong file name'
            assert test_storage.fname_on_disk is None, 'wrong fname on disk'
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
def test_run_by_tap_query_2():
    test_obs_id = 'GS-2017A-Q-58-66-027'
    test_f_id = 'S20170905S0318'
    _write_state()
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        # execution
        with patch('caom2pipe.execute_composable._do_one') as run_mock, \
                patch('gem2caom2.work.TapNoPreviewQuery') as query_mock:
            query_mock.return_value.todo.return_value = [test_obs_id]
            composable.run_by_tap_query()
            assert run_mock.called, 'should have been called'
            args, kwargs = run_mock.call_args
            assert args[3] == 'gem2caom2', 'wrong command'
            test_storage = args[2]
            assert isinstance(test_storage, gem_name.GemName), type(test_storage)
            assert test_storage.obs_id == test_obs_id, 'wrong obs id'
            assert test_storage.file_name is None, 'wrong file name'
            assert test_storage.fname_on_disk is None, 'wrong fname on disk'
            assert test_storage.url is None, 'wrong url'
            assert test_storage.lineage == \
                   '{}/gemini:GEM/{}.fits'.format(test_f_id, test_f_id), \
                'wrong lineage'
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
            assert isinstance(test_storage, gem_name.GemName), type(test_storage)
            assert test_storage.obs_id == test_obs_id, 'wrong obs id'
            assert test_storage.file_name is None, 'wrong file name'
            assert test_storage.fname_on_disk is None, 'wrong fname on disk'
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
@patch('gem2caom2.scrape.read_json_file_list_page')
@patch('gem2caom2.external_metadata.get_obs_metadata')
@patch('caom2pipe.manage_composable.read_obs_from_file')
def test_run_by_builder(read_mock, obs_md_mock, scrape_mock, data_client_mock,
                        repo_mock, exec_mock):
    data_client_mock.return_value.get_file_info.side_effect = \
        gem_mocks.mock_get_file_info
    data_client_mock.return_value.get_file.side_effect = Mock()
    exec_mock.side_effect = Mock()
    repo_mock.return_value.create.side_effect = gem_mocks.mock_repo_create
    repo_mock.return_value.read.side_effect = gem_mocks.mock_repo_read
    repo_mock.return_value.update.side_effect = gem_mocks.mock_repo_update

    def _scrape_mock(ignore1, ignore2):
        return gem_mocks.TEST_TODO_LIST
    scrape_mock.side_effect = _scrape_mock

    obs_md_mock.side_effect = gem_mocks.mock_get_obs_metadata

    def _read_mock(ignore_fqn):
        return SimpleObservation(collection='TEST',
                                 observation_id='TEST_OBS_ID',
                                 algorithm=Algorithm('exposure'))
    read_mock.side_effect = _read_mock

    _write_cert()
    prior_s = datetime.utcnow().timestamp() - 1440 * 60
    _write_state(prior_s)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=f'{gem_mocks.TEST_DATA_DIR}/edu_query')
    try:
        # execution
        test_result = composable._run_by_builder()
        assert test_result == 0, 'wrong result'
    finally:
        os.getcwd = getcwd_orig

    assert repo_mock.return_value.create.called, 'create not called'
    assert repo_mock.return_value.read.called, 'read not called'
    assert exec_mock.called, 'exec mock not called'
    param, level_as = ec.CaomExecute._specify_logging_level_param(logging.ERROR)
    exec_mock.assert_called_with(
        ('gem2caom2 --quiet --cert /usr/src/app/cadcproxy.pem '
         '--in /usr/src/app/logs/GN-2019B-ENG-1-160-008.fits.xml '
         '--out /usr/src/app/logs/GN-2019B-ENG-1-160-008.fits.xml '
         '--external_url '
         'https://archive.gemini.edu/fullheader/N20191101S0007.fits '
         '--plugin '
         '/usr/local/lib/python3.6/site-packages/gem2caom2/gem2caom2.py '
         '--module '
         '/usr/local/lib/python3.6/site-packages/gem2caom2/gem2caom2.py '
         '--lineage N20191101S0007/gemini:GEM/N20191101S0007.fits'),
        level_as), \
        'exec mock wrong parameters'
    assert not data_client_mock.return_value.get_file_info.called, \
        'data client mock get file info should not be not called'
    assert obs_md_mock.called, 'obs mock not called'
    assert scrape_mock.called, 'scrape mock not called'


@patch('sys.exit', Mock(return_value=MyExitError))
@patch('caom2pipe.manage_composable.exec_cmd')
@patch('caom2pipe.execute_composable.CAOM2RepoClient')
@patch('caom2pipe.execute_composable.CadcDataClient')
@patch('caom2pipe.astro_composable.query_tap')
@patch('caom2pipe.manage_composable.query_endpoint')
def test_run_by_public(scrape_mock, tap_mock, data_client_mock,
                       repo_mock, exec_mock):
    data_client_mock.return_value.get_file_info.side_effect = \
        gem_mocks.mock_get_file_info
    data_client_mock.return_value.get_file.side_effect = Mock()
    exec_mock.side_effect = Mock()
    repo_mock.return_value.create.side_effect = gem_mocks.mock_repo_create
    repo_mock.return_value.read.side_effect = gem_mocks.mock_repo_read
    repo_mock.return_value.update.side_effect = gem_mocks.mock_repo_update

    def _tap_mock(ignore1, ignore2):
        c = Column(name='observationID', data=['GN-2019B-ENG-1-160-008'])
        t = Table()
        t.add_column(c)
        return t
    tap_mock.side_effect = _tap_mock

    def _scrape_mock(ignore1):
        result = gem_mocks.Object()

        def _mock_json_return():
            with open(PUBLIC_TEST_JSON) as f:
                stuff = json.loads(f.read())
            return stuff
        result.json = _mock_json_return

        return result
    scrape_mock.side_effect = _scrape_mock

    if not os.path.exists(
            f'/usr/src/app/logs/{gem_mocks.TEST_BUILDER_OBS_ID}.fits.xml'):
        shutil.copy(
            f'{gem_mocks.TEST_DATA_DIR}/expected.xml',
            f'/usr/src/app/logs/{gem_mocks.TEST_BUILDER_OBS_ID}.fits.xml')

    _write_cert()
    prior_s = datetime.utcnow().timestamp() - 1440 * 60
    _write_state(prior_s)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=f'{gem_mocks.TEST_DATA_DIR}/edu_query')
    try:
        # execution
        test_result = composable._run_by_public()
        assert test_result == 0, 'wrong result'
    finally:
        os.getcwd = getcwd_orig

    assert repo_mock.return_value.update.called, 'update not called'
    assert repo_mock.return_value.read.called, 'read not called'
    assert exec_mock.called, 'exec mock not called'
    param, level_as = ec.CaomExecute._specify_logging_level_param(logging.ERROR)
    exec_mock.assert_called_with(
        ('gem2caom2 --quiet --cert /usr/src/app/cadcproxy.pem '
         '--in /usr/src/app/logs/GN-2019B-ENG-1-160-008.fits.xml '
         '--out /usr/src/app/logs/GN-2019B-ENG-1-160-008.fits.xml '
         '--external_url '
         'https://archive.gemini.edu/fullheader/N20191101S0007.fits '
         '--plugin '
         '/usr/local/lib/python3.6/site-packages/gem2caom2/gem2caom2.py '
         '--module '
         '/usr/local/lib/python3.6/site-packages/gem2caom2/gem2caom2.py '
         '--lineage N20191101S0007/gemini:GEM/N20191101S0007.fits'),
        level_as), \
        'exec mock wrong parameters'
    assert not data_client_mock.return_value.get_file_info.called, \
        'data client mock get file info should not be called, file_name is None'
    assert tap_mock.called, 'tap mock not called'
    assert scrape_mock.called, 'scrape mock not called'


def _write_todo(test_obs_id):
    with open(TODO_FILE, 'w') as f:
        f.write('{}\n'.format(test_obs_id))


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
    test_start_time = datetime.fromtimestamp(prior_s)
    test_bookmark = {'bookmarks': {'gemini_timestamp':
                                       {'last_record': test_start_time}}}
    mc.write_as_yaml(test_bookmark, STATE_FILE)


def _write_rejected(test_obs_id):
    content = {'bad_metadata': [test_obs_id]}
    mc.write_as_yaml(content, REJECTED_FILE)


def _write_cert():
    if not os.path.exists('/usr/src/app/cadcproxy.pem'):
        with open('/usr/src/app/cadcproxy.pem', 'w') as f:
            f.write('cadc proxy content')
