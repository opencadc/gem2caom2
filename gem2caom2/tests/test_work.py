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
import sys

from datetime import datetime
from mock import patch, Mock
import pytest

import gem_mocks

from caom2pipe import manage_composable as mc
from gem2caom2 import main_app, composable, gem_name, work

STATE_FILE = '/usr/src/app/state.yml'
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')

multiple_call_count = 0


def test_run_edu_parse():
    start_time = datetime.strptime(
        '2019-10-10T05:04:24.000', mc.ISO_8601_FORMAT)
    end_time = datetime.strptime(
        '2019-10-10T06:04:24.000', mc.ISO_8601_FORMAT)
    # execution
    html_string = open(gem_mocks.FIRST_FILE_LIST, 'r').read()
    work_list_result, max_date_result =\
        work.ArchiveGeminiEduQuery.parse_ssummary_page(
            html_string, start_time, end_time)
    assert work_list_result is not None, 'expected result'
    assert max_date_result is not None, 'expected result'
    assert len(work_list_result) == 1, 'wrong number of results'
    expected_key = 'S20191010S0030'
    assert work_list_result[expected_key].obs_id == \
        'GS-2019B-Q-222-181-001', 'wrong obs id'
    assert work_list_result[expected_key].file_id == expected_key, \
        'wrong file id'


@patch('gem2caom2.composable._get_utcnow')
def test_run_edu_parse_too_many_files(utc_now_patch):
    start_time = datetime.strptime(
        '2003-01-05T00:00:00.000', mc.ISO_8601_FORMAT)
    end_time = datetime.strptime(
        '2003-01-07T00:00:00.000', mc.ISO_8601_FORMAT)
    gem_mocks.mock_write_state(start_time)
    # execution
    with patch('caom2pipe.manage_composable.query_endpoint') as \
            query_endpoint_mock, \
            patch('caom2pipe.execute_composable._storage_name_middle') as \
            run_mock:
        query_endpoint_mock.side_effect = gem_mocks.mock_query_endpoint
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=f'{TEST_DATA_DIR}')
        utc_now_patch.return_value = end_time
        try:
            with pytest.raises(mc.CadcException):
                # raise because 2500 records
                sys.argv = ['test_command']
                composable._run_by_edu_query()
        finally:
            os.getcwd = getcwd_orig

        assert run_mock.called, 'run mock should be called'
        assert run_mock.call_count == 1, 'query limit failure'


class MyTimingException(Exception):
    pass


def test_run_edu_multiple_returns():
    start_time = datetime.strptime(
        '2019-10-10T05:09:24.000', mc.ISO_8601_FORMAT)
    # 24 hours == 1440 minutes
    gem_mocks.mock_write_state(start_time)
    # execution
    with patch('caom2pipe.manage_composable.query_endpoint') as \
            query_endpoint_mock, \
            patch('caom2pipe.execute_composable._do_one') as run_mock:
        query_endpoint_mock.side_effect = _multiple_returns
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=f'{TEST_DATA_DIR}/edu_query')
        try:
            sys.argv = ['test_command']
            composable._run_by_edu_query()
        except MyTimingException as e:
            # how to end the test without trailing timeboxes
            assert run_mock.called, 'should have been called'
            args, kwargs = run_mock.call_args
            assert args[3] == main_app.APPLICATION, 'wrong command'
            test_storage = args[2]
            assert isinstance(
                test_storage, gem_name.GemName), type(test_storage)
            assert run_mock.call_count == 3, 'wrong call count'

        finally:
            os.getcwd = getcwd_orig


def _multiple_returns(url, timeout=-1):
    global multiple_call_count
    result = gem_mocks.Object()
    result.text = ''
    file_pre = url.split('/')[-1].replace('filepre=', '')
    first = f'{TEST_DATA_DIR}/edu_query/first.html'
    second = f'{TEST_DATA_DIR}/edu_query/second.html'
    third = f'{TEST_DATA_DIR}/edu_query/third.html'
    if (multiple_call_count == 0 and file_pre.startswith('S') and
            file_pre.endswith('0') and '20191010' in url):
        multiple_call_count += 1
        with open(first, 'r') as f:
            result.text = f.read()
    elif (multiple_call_count == 1 and file_pre.startswith('N') and
          file_pre.endswith('3') and '20191011' in url):
        multiple_call_count += 1
        with open(second, 'r') as f:
            result.text = f.read()
    elif (multiple_call_count == 2 and file_pre.startswith('S') and
          file_pre.endswith('1') and '20191012' in url):
        multiple_call_count += 1
        with open(third, 'r') as f:
            result.text = f.read()
    else:
        if multiple_call_count > 2:
            import logging
            logging.error(f'wut {url} multiple call count {multiple_call_count}')
            raise MyTimingException(f'url::{url} {multiple_call_count}')
    return result
