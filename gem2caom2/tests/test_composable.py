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

import logging
import os
import pytest
import sys
import traceback

from datetime import datetime
from mock import patch, Mock

from caom2pipe import manage_composable as mc
from gem2caom2 import composable, GemName


PY_VERSION = '3.6'
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
STATE_FILE = '/usr/src/app/state.yml'
TODO_FILE = '/usr/src/app/todo.txt'


class MyExitError(Exception):
    pass


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support single version')
@patch('sys.exit', Mock(return_value=MyExitError))
def test_run_query():
    # preconditions
    now = datetime.utcnow()
    test_start_time = datetime(year=now.year, month=now.month,
                               day=now.day, hour=now.hour,
                               minute=now.minute - 1, second=now.second,
                               microsecond=now.microsecond,
                               tzinfo=now.tzinfo)
    test_bookmark = {'bookmarks': {'gemini_timestamp':
                                   {'last_record': test_start_time}}}
    mc.write_as_yaml(test_bookmark, STATE_FILE)
    start_time = os.path.getmtime(STATE_FILE)

    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=TEST_DATA_DIR)

    try:
        # execution
        with patch('caom2pipe.astro_composable.query_tap') as \
                query_endpoint_mock, \
                patch('caom2pipe.execute_composable.run_by_file_prime') \
                as run_mock:
            query_endpoint_mock.return_value = {
                'observationID': [b'GS-2004A-Q-6-27-0255']}
            composable.run_query()
            assert run_mock.called, 'should have been called'
        end_time = os.path.getmtime(STATE_FILE)
        assert end_time > start_time, 'no execution'
    except Exception as e:
        logging.error(traceback.format_exc())
        assert False
    finally:
        os.getcwd = getcwd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support single version')
@patch('sys.exit', Mock(return_value=MyExitError))
def test_run():
    test_obs_id = 'GS-2004A-Q-6-27-0255'
    test_f_id = '2004may19_0255'
    _write_todo(test_obs_id)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=TEST_DATA_DIR)
    try:
        # execution
        with patch('caom2pipe.execute_composable._do_one') \
                as run_mock:
            composable.run()
            assert run_mock.called, 'should have been called'
            args, kwargs = run_mock.call_args
            assert args[3] == 'gem2caom2', 'wrong command'
            test_storage = args[2]
            assert isinstance(test_storage, GemName), type(test_storage)
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
    except Exception as e:
        logging.error(traceback.format_exc())
        assert False
    finally:
        os.getcwd = getcwd_orig


# gem2caom2 --verbose --cert /usr/src/app/cadcproxy.pem --observation GEMINI
# TX20131117_flt.3002 --out /usr/src/app/logs/TX20131117_flt.3002.fits.xml
# --plugin /usr/local/lib/python3.6/site-packages/gem2caom2/gem2caom2.py
# --module /usr/local/lib/python3.6/site-packages/gem2caom2/gem2caom2.py
# --lineage TX20131117_flt.3002/gemini:GEM/TX20131117_flt.3002.fits
@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support single version')
@patch('sys.exit', Mock(return_value=MyExitError))
def test_run_errors():
    test_obs_id = 'TX20131117_flt.3002'
    test_f_id = 'TX20131117_flt.3002'
    _write_todo(test_obs_id)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=TEST_DATA_DIR)
    try:
        # execution
        with patch('caom2pipe.execute_composable._do_one') \
                as run_mock:
            composable.run()
            assert run_mock.called, 'should have been called'
            args, kwargs = run_mock.call_args
            assert args[3] == 'gem2caom2', 'wrong command'
            test_storage = args[2]
            assert isinstance(test_storage, GemName), type(test_storage)
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
    except Exception as e:
        logging.error(traceback.format_exc())
        assert False
    finally:
        os.getcwd = getcwd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support single version')
@patch('sys.exit', Mock(return_value=MyExitError))
def test_run_query():
    test_obs_id = 'GS-2017A-Q-58-66-027'
    test_f_id = 'S20170905S0318'
    start_time = os.path.getmtime(STATE_FILE)
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=TEST_DATA_DIR)
    try:
        # execution
        with patch('caom2pipe.execute_composable._do_one') as run_mock, \
                patch('gem2caom2.work.read_obs_ids_from_caom') as query_mock:
            query_mock.return_value = [test_obs_id]
            composable.run_query()
            assert run_mock.called, 'should have been called'
            args, kwargs = run_mock.call_args
            assert args[3] == 'gem2caom2', 'wrong command'
            test_storage = args[2]
            assert isinstance(test_storage, GemName), type(test_storage)
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
            test_config = args[0]
            assert isinstance(test_config, mc.Config), 'wrong arg type'
            assert test_config.state_fqn == STATE_FILE, 'wrong state file'
            end_time = os.path.getmtime(STATE_FILE)
            assert end_time > start_time, 'state file not updated'
    except Exception as e:
        logging.error(traceback.format_exc())
        assert False
    finally:
        os.getcwd = getcwd_orig


def _write_todo(test_obs_id):
    with open(TODO_FILE, 'w') as f:
        f.write('{}\n'.format(test_obs_id))
