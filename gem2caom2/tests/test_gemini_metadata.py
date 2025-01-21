# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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

import os

from astropy.io import fits
from astropy.table import Table

from mock import ANY, patch, Mock

from cadcutils import exceptions
from caom2utils.data_util import get_local_file_headers
from caom2pipe import manage_composable as mc
from gem2caom2 import gemini_metadata, gem_name

import gem_mocks


@patch('caom2utils.data_util.get_local_file_headers', autospec=True)
@patch('caom2pipe.client_composable.query_tap_client', autospec=True)
def test_provenance_finder(caom2_mock, local_mock):
    test_file_id = 'rN20123456S9876'
    test_uri = f'gemini:GEMINI/{test_file_id}.fits'
    repaired_data_label = 'GN-2012A-B-012-345-6'
    test_data_label = f'{repaired_data_label}-R'

    def _caom2_mock(ignore1, ignore2):
        return Table.read(
            f'observationID,instrument_name\n{test_data_label},GMOS\n'.split('\n'),
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
    test_config.proxy_fqn = os.path.join(gem_mocks.TEST_DATA_DIR, 'cadcproxy.pem')
    test_config.tap_id = 'ivo://cadc.nrc.ca/test'

    try:
        for test_use_local in [True, False]:
            for test_connected in [True, False]:
                test_config.use_local_files = test_use_local
                if test_connected:
                    test_config.task_types = [mc.TaskType.VISIT]
                else:
                    test_config.task_types = [mc.TaskType.SCRAPE]

                test_subject = gemini_metadata.ProvenanceFinder(Mock(), test_config)
                assert test_subject is not None, (
                    f'ctor does not work:: local {test_use_local}, connected {test_connected}'
                )
                test_result = test_subject.get(test_uri)
                assert test_result is not None, (
                    f'expect a result local {test_use_local}, connected {test_connected}'
                )
                assert test_result == repaired_data_label, (
                    f'data_label should be {repaired_data_label} '
                    f'local {test_use_local}, '
                    f'connected {test_connected}'
                )
    finally:
        os.path.exists = os_path_exists_orig


@patch('caom2pipe.client_composable.ClientCollection')
@patch('gem2caom2.gemini_metadata.retrieve_gemini_headers')
def test_header_not_at_cadc_no_reader(retrieve_gemini_mock, clients_mock, test_config):
    # the file is private, re-ingestion fails to find headers at CADC, needs to go back to archive.gemini.edu
    test_f_name = 'N20220314S0229.fits.bz2'
    test_obs_id = 'GN-CAL20220314-18-083'
    retrieve_gemini_mock.side_effect = gem_mocks._mock_retrieve_headers
    clients_mock.data_client.get_head.side_effect = exceptions.UnexpectedException
    test_storage_name = gem_name.GemName(file_name=test_f_name)
    test_storage_name.obs_id = test_obs_id
    test_result = gemini_metadata.retrieve_headers(test_storage_name, 0, Mock(), clients_mock, test_config)
    assert test_result is not None, 'expect a result'
    assert retrieve_gemini_mock.called, 'retrieve mock not called'
    retrieve_gemini_mock.assert_called_with(test_f_name, ANY, ANY), 'wrong mock args'


@patch('gem2caom2.composable.GemClientCollection')
def test_header_not_at_cadc_no_reader_session_mock(clients_mock, test_data_dir, test_config):
    # the file is private, re-ingestion fails to find headers at CADC, needs to go back to archive.gemini.edu
    test_f_name = 'N20220314S0229.fits.bz2'
    test_obs_id = 'GN-CAL20220314-18-083'
    clients_mock.data_client.get_head.side_effect = exceptions.UnexpectedException
    session_return = gem_mocks.Object()
    with open (f'{test_data_dir}/N20220314S0229.fits') as f_in:
        session_return.text = f_in.read()
    clients_mock.gemini_session.get.return_value = session_return
    test_storage_name = gem_name.GemName(file_name=test_f_name)
    test_storage_name.obs_id = test_obs_id
    test_result = gemini_metadata.retrieve_headers(test_storage_name, 0, Mock(), clients_mock, test_config)
    assert test_result is not None, 'expect a result'
    assert clients_mock.gemini_session.get.called, 'get mock not called'
    clients_mock.gemini_session.get.assert_called_with(
        'https://archive.gemini.edu/fullheader/N20220314S0229.fits.bz2', timeout=20
    ), 'wrong mock args'
