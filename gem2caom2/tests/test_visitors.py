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

import os
import pytest
import shutil

from datetime import datetime, timezone
from mock import patch, Mock

from cadcutils import exceptions
from cadcdata import FileInfo
from caom2 import Instrument
from gem2caom2 import ghost_preview_augmentation, preview_augmentation, pull_augmentation, cleanup_augmentation
from gem2caom2 import gemini_metadata, util
from gem2caom2 import svofps, gem_name
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
import gem_mocks

pytest.main(args=['-s', os.path.abspath(__file__)])
TEST_OBS = 'GN-2013B-Q-28-150-002'
TEST_FILE = 'N20131203S0006.jpg'
TEST_FP_OBS = 'GN-2015A-C-1-20-001'
TEST_FP_FILE = 'N20150404S0726.fits'
TEST_PRODUCT_ID = 'GN2001BQ013-04'


def test_preview_augment_known_no_preview(test_data_dir, test_config, tmp_path):
    # rejected file exists that says there's a preview known to not
    # exist, so trying to generate a thumbnail will result in no
    # change to the plane/artifact structure

    test_config.change_working_directory(tmp_path.as_posix())
    obs = mc.read_obs_from_file(f'{test_data_dir}/visit_obs_start.xml')
    obs.planes[TEST_PRODUCT_ID].data_release = datetime.now(tz=timezone.utc).replace(tzinfo=None)
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'

    if os.path.exists(test_config.rejected_fqn):
        os.unlink(test_config.rejected_fqn)
    test_rejected = mc.Rejected(test_config.rejected_fqn)
    test_rejected.record(mc.Rejected.NO_PREVIEW, f'{TEST_PRODUCT_ID}.jpg')
    test_reporter = mc.ExecutionReporter2(test_config)
    test_reporter._observable._rejected = test_rejected
    test_storage_name = gem_name.GemName(file_name=TEST_FP_FILE)
    cadc_client_mock = Mock()
    clients_mock = Mock()
    clients_mock.data_client = cadc_client_mock
    kwargs = {
        'working_directory': test_data_dir,
        'clients': clients_mock,
        'stream': 'stream',
        'reporter': test_reporter,
        'storage_name': test_storage_name,
    }

    with patch('caom2pipe.manage_composable.http_get') as http_mock, patch(
        'caom2pipe.manage_composable.get_artifact_metadata'
    ) as art_mock, patch('caom2pipe.manage_composable.exec_cmd') as exec_mock:
        cadc_client_mock.return_value.data_get.return_value = mc.CadcException('test')
        obs = preview_augmentation.visit(obs, **kwargs)
        assert not http_mock.called, 'http mock should not be called'
        assert not cadc_client_mock.put.called, 'put mock should not be called'
        assert not art_mock.called, 'art mock should not be called'
        assert not exec_mock.called, 'exec mock should not be called'
        assert obs is not None, 'expect a result'
        assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'no new artifacts'

    test_rejected.persist_state()
    assert os.path.exists(test_config.rejected_fqn)


def test_preview_augment_unknown_no_preview(test_data_dir, test_config, tmp_path):
    # what happens when it's not known that there's no preview
    test_config.change_working_directory(tmp_path)
    obs = mc.read_obs_from_file(f'{test_data_dir}/visit_obs_start.xml')
    obs.planes[TEST_PRODUCT_ID].data_release = datetime.now(tz=timezone.utc).replace(tzinfo=None)
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'

    test_reporter = mc.ExecutionReporter2(test_config)
    test_storage_name = gem_name.GemName(file_name=f'{TEST_PRODUCT_ID}.fits')

    cadc_client_mock = Mock()
    clients_mock = Mock()
    clients_mock.data_client = cadc_client_mock
    kwargs = {
        'working_directory': test_data_dir,
        'clients': clients_mock,
        'stream': 'stream',
        'reporter': test_reporter,
        'storage_name': test_storage_name,
    }

    with patch(
        'caom2pipe.manage_composable.http_get',
        side_effect=mc.CadcException('Not Found for url: https://archive.gemini.edu/preview'),
    ) as http_mock, patch('caom2pipe.manage_composable.get_artifact_metadata') as art_mock, patch(
        'caom2pipe.manage_composable.exec_cmd'
    ) as exec_mock:
        cadc_client_mock.get.side_effect = exceptions.UnexpectedException('test')
        obs = preview_augmentation.visit(obs, **kwargs)
        assert obs is not None, 'expect result'
        test_url = f'{preview_augmentation.PREVIEW_URL}{TEST_PRODUCT_ID}.fits'
        test_prev = f'{test_data_dir}/{TEST_PRODUCT_ID}.jpg'
        http_mock.assert_called_with(test_url, test_prev), 'mock not called'
        assert not cadc_client_mock.put.called, 'put mock should not be called'
        assert not art_mock.called, 'art mock should not be called'
        assert not exec_mock.called, 'exec mock should not be called'


@patch('caom2utils.data_util.get_file_type')
@patch('gem2caom2.gemini_metadata.retrieve_json')
@patch('caom2pipe.manage_composable.http_get')
def test_pull_augmentation(http_mock, json_mock, file_type_mock, test_config, test_data_dir, tmp_path):
    test_config.change_working_directory(tmp_path)
    obs = mc.read_obs_from_file(f'{test_data_dir}/visit_original_uri_start.xml')
    obs.planes[TEST_PRODUCT_ID].data_release = datetime.now(tz=timezone.utc).replace(tzinfo=None)
    original_uri = 'gemini:GEMINI/GN2001BQ013-04.fits'
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'
    assert original_uri in obs.planes[TEST_PRODUCT_ID].artifacts.keys(), 'initial condition'
    test_uri = f'{test_config.scheme}:{test_config.collection}/{TEST_PRODUCT_ID}.fits'

    test_reporter = mc.ExecutionReporter2(test_config)
    cadc_client_mock = Mock()
    cadc_client_mock.return_value.info.return_value = None
    clients_mock = Mock()
    clients_mock.data_client = cadc_client_mock
    json_mock.side_effect = gem_mocks.mock_retrieve_json
    filter_cache = svofps.FilterMetadataCache(Mock())
    test_storage_name = gem_name.GemName(file_name='GN2001BQ013-04.fits', filter_cache=filter_cache)
    file_type_mock.return_values = 'application/fits'
    kwargs = {
        'working_directory': test_data_dir,
        'clients': clients_mock,
        'reporter': test_reporter,
        'storage_name': test_storage_name,
    }

    test_precondition = gemini_metadata.GeminiMetaVisitRunnerMeta(clients_mock, test_config, [], test_reporter)
    test_precondition._storage_name = test_storage_name
    test_precondition._set_preconditions()

    obs = pull_augmentation.visit(obs, **kwargs)
    test_url = f'{pull_augmentation.FILE_URL}/{TEST_PRODUCT_ID}.fits'
    test_prev = f'{test_data_dir}/{TEST_PRODUCT_ID}.fits'
    http_mock.assert_called_with(test_url, test_prev), 'mock not called'
    assert cadc_client_mock.put.called, 'put mock not called'
    cadc_client_mock.put.assert_called_with(test_data_dir, 'gemini:GEMINI/GN2001BQ013-04.fits'), 'wrong put args'
    assert obs is not None, 'expect a result'
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'no new artifacts'
    try:
        ignore = obs.planes[TEST_PRODUCT_ID].artifacts[test_uri]
    except KeyError as ke:
        # because CAOM does magic
        result = obs.planes[TEST_PRODUCT_ID].artifacts[original_uri]
        assert result.uri == test_uri, f'wrong uri {result.uri}'


def test_preview_augment_delete_preview(test_data_dir, test_config, tmp_path):
    test_config.change_working_directory(tmp_path.as_posix())
    # plane starts with a preview artifact, but it represents a non-existent
    # file, so remove the artifact from the CAOM observation
    test_product_id = 'S20080610S0045'
    fqn = os.path.join(test_data_dir, 'GS-2008A-C-5-35-002.fits.xml')
    obs = mc.read_obs_from_file(fqn)
    assert len(obs.planes[test_product_id].artifacts) == 2, 'initial condition'
    test_rejected = mc.Rejected('/tmp/nonexistent')
    test_rejected.content = {
        'bad_metadata': [],
        'no_preview': [
            'S20080610S0043.jpg',
            'S20080610S0041.jpg',
            'S20080610S0044.jpg',
            'S20080610S0045.jpg',
        ],
    }
    test_reporter = mc.ExecutionReporter2(test_config)
    test_reporter._observable._rejected = test_rejected
    test_storage_name = gem_name.GemName(file_name=f'{test_product_id}.fits')
    kwargs = {
        'working_directory': test_data_dir,
        'clients': None,
        'stream': 'stream',
        'reporter': test_reporter,
        'storage_name': test_storage_name,
    }
    obs = preview_augmentation.visit(obs, **kwargs)
    assert obs is not None, 'expect a result'
    assert len(obs.planes[test_product_id].artifacts) == 1, 'post condition'


@patch('caom2pipe.manage_composable.http_get')
def test_preview_augment(http_mock, test_data_dir, test_config, tmp_path):
    test_config.change_working_directory(tmp_path.as_posix())
    # this should result in two new artifacts being added to the plane
    # one for a thumbnail and one for a preview

    obs = mc.read_obs_from_file(f'{test_data_dir}/visit_obs_start.xml')
    obs.planes[TEST_PRODUCT_ID].data_release = datetime.now(tz=timezone.utc).replace(tzinfo=None)
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'

    test_reporter = mc.ExecutionReporter2(test_config)
    cadc_client_mock = Mock()
    clients_mock = Mock()
    clients_mock.data_client = cadc_client_mock
    test_storage_name = gem_name.GemName(file_name=f'{TEST_PRODUCT_ID}.fits')
    kwargs = {
        'working_directory': '/test_files',
        'clients': clients_mock,
        'reporter': test_reporter,
        'storage_name': test_storage_name,
    }

    test_prev = f'/test_files/{TEST_PRODUCT_ID}.jpg'
    if os.path.exists(test_prev):
        os.unlink(test_prev)

    def _get_mock(url_ignore, fqn_ignore):
        shutil.copy(f'{test_data_dir}/{TEST_FILE}', f'/test_files/{TEST_PRODUCT_ID}.jpg')

    try:
        cadc_client_mock.get.side_effect = exceptions.UnexpectedException('test')
        http_mock.side_effect = _get_mock
        obs = preview_augmentation.visit(obs, **kwargs)
        test_url = f'{preview_augmentation.PREVIEW_URL}{TEST_PRODUCT_ID}.fits'
        assert http_mock.called, 'http mock should be called'
        http_mock.assert_called_with(test_url, test_prev), 'mock not called'
        assert cadc_client_mock.put.called, 'put mock not called'
        cadc_client_mock.put.assert_called_with(
            '/test_files',
            'cadc:GEMINI/GN2001BQ013-04_th.jpg',
        ), 'wrong put arguments'
        assert obs is not None, 'expect a result'
        assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 3, 'two new artifacts'
        prev_uri = mc.build_uri(test_config.collection, f'{TEST_PRODUCT_ID}.jpg', test_config.scheme)
        thumb_uri = mc.build_uri(test_config.collection, f'{TEST_PRODUCT_ID}_th.jpg', test_config.preview_scheme)
        assert prev_uri in obs.planes[TEST_PRODUCT_ID].artifacts.keys(), 'no preview'
        assert thumb_uri in obs.planes[TEST_PRODUCT_ID].artifacts, 'no thumbnail'
    finally:
        if os.path.exists(test_prev):
            os.unlink(test_prev)


@patch('caom2pipe.manage_composable.http_get')
def test_preview_augment_failure(http_mock, test_data_dir, test_config, tmp_path):
    test_config.change_working_directory(tmp_path.as_posix())
    # mimic 'Not Found' behaviour
    # this should result in no new artifacts being added to the plane
    # but a record for 'no preview exists at Gemini' added to the
    # record

    def _failure_mock(ignore_url, ignore_local_fqn):
        raise mc.CadcException(
            'Could not retrieve /usr/src/app/N20211007A0003/'
            'N20211007A0003b.jpg from '
            'https://archive.gemini.edu/preview/N20211007A0003b.fits. Failed '
            'with 404 Client Error: Not Found for url: '
            'https://archive.gemini.edu/preview/N20211007A0003b.fits'
        )

    obs = mc.read_obs_from_file(f'{test_data_dir}/visit_obs_start.xml')
    obs.planes[TEST_PRODUCT_ID].data_release = datetime.now(tz=timezone.utc).replace(tzinfo=None)
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'

    test_reporter = mc.ExecutionReporter2(test_config)
    cadc_client_mock = Mock()
    clients_mock = Mock()
    clients_mock.data_client = cadc_client_mock
    test_storage_name = gem_name.GemName(file_name=f'{TEST_PRODUCT_ID}.fits')
    kwargs = {
        'working_directory': '/test_files',
        'clients': clients_mock,
        'reporter': test_reporter,
        'storage_name': test_storage_name,
    }

    test_prev = f'/test_files/{TEST_PRODUCT_ID}.jpg'
    if os.path.exists(test_prev):
        os.unlink(test_prev)

    try:
        cadc_client_mock.get.side_effect = exceptions.UnexpectedException('test')
        http_mock.side_effect = _failure_mock
        obs = preview_augmentation.visit(obs, **kwargs)
        test_url = f'{preview_augmentation.PREVIEW_URL}{TEST_PRODUCT_ID}.fits'
        assert http_mock.called, 'http mock should be called'
        http_mock.assert_called_with(test_url, test_prev), 'mock not called'
        assert not cadc_client_mock.put.called, 'put mock should not be called'
        assert obs is not None, 'expect a result'
        assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'same as the pre-condition'
        prev_uri = mc.build_uri(test_config.collection, f'{TEST_PRODUCT_ID}.jpg', test_config.scheme)
        thumb_uri = mc.build_uri(test_config.collection, f'{TEST_PRODUCT_ID}_th.jpg', test_config.preview_scheme)
        assert prev_uri not in obs.planes[TEST_PRODUCT_ID].artifacts.keys(), 'should be no preview'
        assert thumb_uri not in obs.planes[TEST_PRODUCT_ID].artifacts, 'should be no thumbnail'
        assert not (test_reporter._observable._rejected.is_no_preview(prev_uri)), 'preview should be tracked'

        assert http_mock.call_count == 1, 'wrong number of calls'
        # now try again to generate the preview, and ensure that the
        # rejected tracking is working
        obs = preview_augmentation.visit(obs, **kwargs)
        assert obs is not None, 'expect a result the second time'
        assert http_mock.call_count == 1, 'never even tried to retrieve it'
    finally:
        if os.path.exists(test_prev):
            os.unlink(test_prev)


def test_cleanup(test_data_dir):
    # test that cleanup occurs where it should
    test_kwargs = {}
    test_cleanup_file = f'{test_data_dir}/cleanup_aug_start_obs.xml'
    test_artifact_id = 'cadc:GEMINI/N20140811S0033_BIAS_th.jpg'
    obs = mc.read_obs_from_file(test_cleanup_file)
    initial_all_artifact_keys = cc.get_all_artifact_keys(obs)
    assert test_artifact_id in initial_all_artifact_keys, 'wrong initial conditions'
    obs = cleanup_augmentation.visit(obs, **test_kwargs)
    assert obs is not None, 'expect a result'
    post_all_artifact_keys = cc.get_all_artifact_keys(obs)
    assert test_artifact_id not in post_all_artifact_keys, 'wrong post conditions'

    # test that cleaning up a clean observation won't break that
    # observation
    obs = cleanup_augmentation.visit(obs, **test_kwargs)

    # test that cleanup doesn't occur where it shouldn't
    test_no_cleanup_file = f'{test_data_dir}/cleanup_no_cleanup_aug_start_obs.xml'
    no_cleanup_obs = mc.read_obs_from_file(test_no_cleanup_file)
    all_artifact_keys = cc.get_all_artifact_keys(no_cleanup_obs)
    assert len(all_artifact_keys) == 6, 'wrong no cleanup initial conditions'
    no_cleanup_obs = cleanup_augmentation.visit(no_cleanup_obs, **test_kwargs)
    all_artifact_keys = cc.get_all_artifact_keys(no_cleanup_obs)
    assert len(all_artifact_keys) == 6, 'wrong no cleanup post conditions, should not be different'

    # test a FOX observation, because it's odd to begin with
    test_fox_file = f'{test_data_dir}/cleanup_fox_aug_start.xml'
    fox_obs = mc.read_obs_from_file(test_fox_file)
    all_artifact_keys = cc.get_all_artifact_keys(fox_obs)
    initial_fox_length = len(all_artifact_keys)
    fox_obs = cleanup_augmentation.visit(fox_obs, **test_kwargs)
    post_fox_length = len(cc.get_all_artifact_keys(fox_obs))
    assert initial_fox_length == post_fox_length, 'wrong fox post conditions'


@patch('caom2pipe.client_composable.ClientCollection')
@patch('caom2pipe.manage_composable.http_get')
def test_look_pull_and_put(http_mock, mock_client, test_data_dir):
    test_storage_name = 'cadc:GEMINI/TEST.fits'
    mock_client.info.return_value = FileInfo(
        id=test_storage_name,
        size=1234,
        md5sum='9473fdd0d880a43c21b7778d34872157',
    )
    f_name = 'TEST.fits'
    url = f'https://localhost/{f_name}'
    test_config = mc.Config()
    test_config.observe_execution = True
    mock_client.info.return_value = None
    test_fqn = os.path.join(test_data_dir, f_name)
    pull_augmentation.look_pull_and_put(
        test_storage_name,
        test_fqn,
        url,
        mock_client,
        'md5:01234',
    )
    mock_client.data_client.put.assert_called_with(test_data_dir, test_storage_name), 'mock not called'
    http_mock.assert_called_with(url, test_fqn), 'http mock not called'


def test_ghost_preview_augmentation(test_config, test_data_dir, tmp_path):
    test_config.change_working_directory(tmp_path.as_posix())
    test_f_id = 'S20230518S0121'
    obs = mc.read_obs_from_file(f'{test_data_dir}/GHOST/{test_f_id}.expected.xml')
    test_reporter = mc.ExecutionReporter2(test_config)
    test_storage_name = gem_name.GemName(file_name=f'/test_files/{test_f_id}.fits')
    test_storage_name._obs_id = 'GS-2023A-SV-101-13-009'
    test_storage_name._destination_uris = [f'gemini:GEMINI/{test_f_id}.fits']
    test_storage_name._product_id = test_f_id
    test_storage_name._file_name = f'{test_f_id}.fits'
    test_storage_name._source_names = [f'/test_files/{test_f_id}.fits']
    kwargs = {
        'working_directory': test_data_dir,
        'clients': None,
        'reporter': test_reporter,
        'storage_name': test_storage_name,
    }
    assert len(obs.planes[test_f_id].artifacts) == 1, 'pre-condition'
    obs.instrument = Instrument(name=util.Inst.GMOS.value)
    obs = ghost_preview_augmentation.visit(obs, **kwargs)
    assert obs is not None, 'expect a result'
    assert len(obs.planes[test_f_id].artifacts) == 1, 'GMOS post-condition'
    obs.instrument = Instrument(name=util.Inst.GHOST.value)
    obs = ghost_preview_augmentation.visit(obs, **kwargs)
    assert obs is not None, 'expect a result'
    assert len(obs.planes[test_f_id].artifacts) == 3, 'GHOST post-condition'


def test_ghost_preview_augmentation_2(test_config, test_data_dir, tmp_path):
    test_config.change_working_directory(tmp_path.as_posix())
    test_f_id = 'S20231208S0060'
    obs = mc.read_obs_from_file(f'{test_data_dir}/GHOST/{test_f_id}.expected.xml')
    test_reporter = mc.ExecutionReporter2(test_config)
    test_storage_name = gem_name.GemName(file_name=f'/test_files/{test_f_id}.fits')
    test_storage_name._obs_id = 'GS-2023B-FT-104-11-001'
    test_storage_name._destination_uris = [f'gemini:GEMINI/{test_f_id}.fits']
    test_storage_name._product_id = test_f_id
    test_storage_name._file_name = f'{test_f_id}.fits'
    test_storage_name._source_names = [f'/test_files/{test_f_id}.fits']
    kwargs = {
        'working_directory': test_data_dir,
        'clients': None,
        'reporter': test_reporter,
        'storage_name': test_storage_name,
    }
    assert len(obs.planes[test_f_id].artifacts) == 1, 'pre-condition'
    obs.instrument = Instrument(name=util.Inst.GHOST.value)
    obs = ghost_preview_augmentation.visit(obs, **kwargs)
    assert obs is not None, 'expect a result'
    assert len(obs.planes[test_f_id].artifacts) == 1, 'GHOST post-condition, there is no preview data'
