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
#  $Revision: 4 $
#
# ***********************************************************************
#

import pytest
from unittest.mock import ANY, Mock, patch
from datetime import datetime, timedelta, timezone
from gem2caom2.gem_name import GemName
from gem2caom2.pull_augmentation import FILE_URL, PullVisitor
from shutil import copyfile


@pytest.fixture
def test_clients():
    clients = Mock()
    clients.data_client.info.return_value = None
    return clients


@pytest.fixture
def test_reporter():
    reporter = Mock()
    reporter._observable.rejected.is_bad_metadata.return_value = False
    return reporter


@pytest.fixture
def test_storage_name():
    storage_name = Mock()
    storage_name.file_name = 'test_uri.fits'
    storage_name.file_uri = 'test_uri.fits'
    storage_name.file_info = {'test_uri.fits': Mock(md5sum='test_md5sum')}
    storage_name.metadata = {}
    return storage_name


@patch('gem2caom2.pull_augmentation.retrieve_headers')
@patch('gem2caom2.pull_augmentation.http_get')
def test_pull_visitor_success(
    http_get_mock,
    retrieve_headers_mock,
    test_clients,
    test_reporter,
    test_config,
    tmp_path,
):
    # file is not at CADC
    test_config.change_working_directory(tmp_path.as_posix())

    test_storage_name = GemName(file_name='S20050825S0143.fits')
    test_storage_name._json_metadata[test_storage_name.file_uri] = {
        'data_md5': 'test_md5sum',
        'release': datetime.now(tz=timezone.utc).replace(tzinfo=None) - timedelta(days=1),
    }
    test_storage_name._file_info = {}

    http_get_mock.side_effect = copyfile(
        '/test_files/S20050825S0143.fits.bz2', f'{tmp_path}/S20050825S0143.fits.bz2'
    )

    visitor = PullVisitor(
        working_directory=tmp_path,
        clients=test_clients,
        reporter=test_reporter,
        storage_name=test_storage_name,
        config=test_config
    )
    visitor.visit()

    http_get_mock.assert_called_once()
    http_get_mock.assert_called_with(
        f'{FILE_URL}/S20050825S0143.fits.bz2', f'{tmp_path}/S20050825S0143.fits.bz2'
    )
    test_clients.data_client.put.assert_called_once()
    test_clients.data_client.put.assert_called_with(tmp_path.as_posix(), test_storage_name.file_uri)
    retrieve_headers_mock.assert_not_called()
    test_clients.data_client.info.assert_called_once()
    test_clients.data_client.info.assert_called_with(test_storage_name.file_uri)


@patch('gem2caom2.pull_augmentation.retrieve_headers')
def test_pull_visitor_bad_metadata(
    retrieve_headers_mock, test_clients, test_reporter, test_storage_name, test_config, tmp_path
):
    test_reporter._observable.rejected.is_bad_metadata.return_value = True

    visitor = PullVisitor(
        working_directory=tmp_path,
        clients=test_clients,
        reporter=test_reporter,
        storage_name=test_storage_name,
        config=test_config
    )
    visitor.visit()

    retrieve_headers_mock.assert_not_called()
    test_clients.data_client.put.assert_not_called()
    test_reporter._observable.rejected.is_bad_metadata.assert_called_with(test_storage_name.obs_id)


@patch('gem2caom2.pull_augmentation.retrieve_headers')
def test_pull_visitor_proprietary_file(retrieve_headers_mock, test_clients, test_reporter, test_config, tmp_path):
    test_config.change_working_directory(tmp_path.as_posix())
    test_storage_name = GemName(file_name='S20050825S0143.fits')
    test_storage_name._json_metadata[test_storage_name.file_uri] = {
        'data_md5': 'test_md5sum',
        'release': datetime.now(tz=timezone.utc).replace(tzinfo=None) + timedelta(days=1),
    }
    test_storage_name._file_info = {}

    visitor = PullVisitor(
        working_directory=tmp_path,
        clients=test_clients,
        reporter=test_reporter,
        storage_name=test_storage_name,
        config=test_config
    )
    visitor.visit()

    retrieve_headers_mock.assert_called_once()
    retrieve_headers_mock.assert_called_with(test_storage_name, 0, ANY, ANY, test_config)
    test_clients.data_client.put.assert_not_called()


@patch('gem2caom2.gemini_metadata.retrieve_gemini_headers')
@patch('gem2caom2.pull_augmentation.http_get')
def test_pull_visitor_file_already_exists(
    http_get_mock, retrieve_headers_mock, test_clients, test_reporter, test_config, tmp_path
):
    # file is at CADC, should be no calls to archive.gemini.edu
    test_config.use_local_files = False
    test_storage_name = GemName(file_name='S20050825S0143.fits')
    test_storage_name._json_metadata[test_storage_name.file_uri] = {
        'data_md5': 'test_md5sum',
        'release': datetime.now(tz=timezone.utc).replace(tzinfo=None) - timedelta(days=1),
    }
    test_storage_name.file_info = {test_storage_name.file_uri: Mock(md5sum='md5:test_md5sum')}
    test_clients.data_client.info.return_value = Mock(md5sum='md5:test_md5sum')

    visitor = PullVisitor(
        working_directory=tmp_path,
        clients=test_clients,
        reporter=test_reporter,
        storage_name=test_storage_name,
        config=test_config
    )
    visitor.visit()

    http_get_mock.assert_not_called()
    retrieve_headers_mock.assert_not_called()
    test_clients.data_client.put.assert_not_called()
    test_clients.data_client.get_head.assert_called()
    test_clients.data_client.get_head.assert_called_with(test_storage_name.file_uri)
