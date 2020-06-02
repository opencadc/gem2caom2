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
import os

from datetime import datetime
from mock import patch, Mock

import gem_mocks

from caom2utils.fits2caom2 import ObsBlueprint
from caom2pipe import manage_composable as mc
from gem2caom2 import work, external_metadata, gem_name, main_app

DATA_FILE = f'{gem_mocks.TEST_DATA_DIR}/edu_query/incremental.json'


def test_file_listing_parse():
    end_time = datetime.utcnow().timestamp()
    start_time = datetime.fromtimestamp(end_time - 2 * 1440 * 60)
    # execution
    with patch('gem2caom2.scrape.read_json_file_list_page') as scrape_mock:
        def _scrape_mock(ignore1, ignore2):
            return gem_mocks.TEST_TODO_LIST
        scrape_mock.side_effect = _scrape_mock

        query_subject = work.FileListingQuery(start_time)
        work_list_result = query_subject.todo()
        assert work_list_result is not None, 'expected result'
        assert len(work_list_result) == 6, 'wrong number of results'
        expected_key = 'N20191101S0002.fits'
        first_result = work_list_result.popitem(last=False)
        assert first_result[1] == expected_key, 'wrong file name'
        assert first_result[0] == 1572566500.951468, 'wrong timestamp'


def test_gemini_incremental_query():
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        test_config = mc.Config()
        test_config.get_executors()
        external_metadata.init_global(incremental=True, config=test_config)
        with patch('caom2pipe.manage_composable.query_endpoint') as query_mock:
            query_mock.side_effect = _mock_endpoint
            test_config = mc.Config()
            test_start_time = datetime(year=2019, month=12, day=1)
            test_end_time = datetime(year=2020, month=12, day=1)
            test_subject = work.GeminiIncrementalQuery(
                datetime.utcnow(), test_config)
            test_result = test_subject.todo(test_start_time, test_end_time)
            assert test_result is not None, 'expect a result'
            assert len(test_result) == 500, 'wrong number of results'
            assert isinstance(test_result[0], gem_name.GemNameBuilder), \
                'wrong result content'
            assert test_result[0].obs_id == 'GN-2001A-C-24-1-073', \
                'wrong obs id'
            assert test_result[0].file_id == '01APR19_073', 'wrong file id'
            assert test_result[0].last_modified_s == 1576541884.426377, \
                'wrong last modified'

            test_blueprint = ObsBlueprint()
            test_uri = mc.build_uri(
                gem_name.ARCHIVE, f'{test_result[0].file_id}.fits.gz')
            main_app.accumulate_fits_bp(
                test_blueprint, test_result[0].file_id, test_uri)
            assert test_blueprint._get('Chunk.time.exposure') == \
                   'get_exposure(header)', 'wrong blueprint'
            assert main_app.get_proposal_id(None) == 'GN-2001A-C-24', \
                'wrong proposal id'
    finally:
        os.getcwd = getcwd_orig


def _mock_endpoint(url, timeout=-1):
    def x():
        with open(DATA_FILE, 'r') as f:
            temp = f.read()
        return json.loads(temp)

    result = gem_mocks.Object()
    result.json = x

    if url.startswith('http://arcdev'):
        assert '2019-12-01T00:00:00.000000%202020-12-01T00:00:00.000000' in \
               url, 'wrong url format'
    else:
        raise mc.CadcException('wut?')
    return result
