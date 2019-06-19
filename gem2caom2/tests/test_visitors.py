# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
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

from datetime import datetime
from mock import patch, Mock

from caom2 import ChecksumURI, Dimension2D, Artifact, ReleaseType, ProductType
from gem2caom2 import preview_augmentation, GemName, SCHEME, ARCHIVE
from caom2pipe import manage_composable as mc

pytest.main(args=['-s', os.path.abspath(__file__)])
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
TEST_OBS = 'GN-2013B-Q-28-150-002'
TEST_FILE = 'N20131203S0006.jpg'
# TEST_FP_OBS = 'GN-2015A-Q-91-5-002'
TEST_FP_OBS = 'GN-2015A-C-1-20-001'
# TEST_FP_FILE = 'N20150216S0142.fits'
TEST_FP_FILE = 'N20150404S0726.fits'
TEST_OBS_FILE = '{}/{}'.format(TEST_DATA_DIR, 'visit_obs_start.xml')
TEST_PRODUCT_ID = 'GN2001BQ013-04'
REJECTED_FILE = os.path.join(TEST_DATA_DIR, 'rejected.yml')


@pytest.mark.skip(reason='Decide what to do about GemName.obs_id value first')
def test_preview_augment_plane(mock_obs_id):
    mock_obs_id.return_value = TEST_OBS
    thumb = os.path.join(TEST_DATA_DIR,
                         'GMOS/{}'.format(GemName(TEST_FILE).thumb))
    if os.path.exists(thumb):
        os.remove(thumb)
    test_fqn = os.path.join(TEST_DATA_DIR, '{}/{}.in.xml'.format(
        'GMOS', TEST_OBS))
    test_obs = mc.read_obs_from_file(test_fqn)
    assert len(test_obs.planes[TEST_OBS].artifacts) == 2
    thumba = '{}:{}/N20131203S0006_th.jpg'.format(SCHEME, ARCHIVE)

    test_kwargs = {'working_directory': '{}/GMOS'.format(TEST_DATA_DIR),
                   'cadc_client': None,
                   'stream': 'default'}
    test_result = preview_augmentation.visit(test_obs, **test_kwargs)
    assert test_result is not None, 'expected a visit return value'
    assert test_result['artifacts'] == 1
    assert len(test_obs.planes[TEST_OBS].artifacts) == 3
    assert os.path.exists(thumb)
    assert test_obs.planes[TEST_OBS].artifacts[thumba].content_checksum == \
        ChecksumURI('md5:a8c106c04db4c148695787bfc364cbd8'), \
        'thumb checksum failure'

    # now do updates
    test_obs.planes[TEST_OBS].artifacts[thumba].content_checksum = \
        ChecksumURI('19661c3c2508ecc22425ee2a05881ed4')
    test_result = preview_augmentation.visit(test_obs, **test_kwargs)
    assert test_result is not None, 'expected update visit return value'
    assert test_result['artifacts'] == 1
    assert len(test_obs.planes) == 1
    assert len(test_obs.planes[TEST_OBS].artifacts) == 3
    assert os.path.exists(thumb)
    assert test_obs.planes[TEST_OBS].artifacts[thumba].content_checksum == \
        ChecksumURI('md5:a8c106c04db4c148695787bfc364cbd8'), \
        'thumb update failed'


def test_preview_augment():
    obs = mc.read_obs_from_file(TEST_OBS_FILE)
    obs.planes[TEST_PRODUCT_ID].data_release = datetime.utcnow()
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'

    test_rejected = mc.Rejected(REJECTED_FILE)
    cadc_client_mock = Mock()
    kwargs = {'working_directory': TEST_DATA_DIR,
              'cadc_client': cadc_client_mock,
              'stream': 'stream',
              'rejected': test_rejected}

    with patch('caom2pipe.manage_composable.http_get') as http_mock, \
            patch('caom2pipe.manage_composable.data_put') as ad_put_mock, \
            patch('caom2pipe.manage_composable.get_artifact_metadata') as \
            art_mock, \
            patch('caom2pipe.manage_composable.exec_cmd') as exec_mock:
        def _art_mock(fqn, product_type, release_type, uri, temp):
            if '_th' in uri:
                return Artifact(uri,
                                ProductType.PREVIEW,
                                ReleaseType.DATA,
                                'image/jpeg',
                                13,
                                ChecksumURI('md5:13'))
            else:
                return Artifact(uri,
                                ProductType.PREVIEW,
                                ReleaseType.DATA,
                                'image/jpeg',
                                12,
                                ChecksumURI('md5:12'))

        cadc_client_mock.return_value.data_get.return_value = mc.CadcException(
            'test')
        art_mock.side_effect = _art_mock
        result = preview_augmentation.visit(obs, **kwargs)
        test_url = '{}{}.fits'.format(preview_augmentation.PREVIEW_URL,
                                      TEST_PRODUCT_ID)
        test_prev = '{}/{}.jpg'.format(TEST_DATA_DIR, TEST_PRODUCT_ID)
        http_mock.assert_called_with(test_url, test_prev),  'mock not called'
        assert ad_put_mock.called, 'ad put mock not called'
        assert art_mock.called, 'art mock not called'
        assert exec_mock.called, 'exec mock not called'
        assert result is not None, 'expect a result'
        assert result['artifacts'] == 2, 'no artifacts should be updated'
        assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 3, \
            'two new artifacts'


def test_preview_augment_known_no_preview():
    try:
        obs = mc.read_obs_from_file(TEST_OBS_FILE)
        obs.planes[TEST_PRODUCT_ID].data_release = datetime.utcnow()
        assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'

        if os.path.exists(REJECTED_FILE):
            os.unlink(REJECTED_FILE)
        test_rejected = mc.Rejected(REJECTED_FILE)
        test_rejected.record(
            mc.Rejected.NO_PREVIEW, '{}.jpg'.format(TEST_PRODUCT_ID))

        cadc_client_mock = Mock()
        kwargs = {'working_directory': TEST_DATA_DIR,
                  'cadc_client': cadc_client_mock,
                  'stream': 'stream',
                  'rejected': test_rejected}

        with patch('caom2pipe.manage_composable.http_get') as http_mock, \
                patch('caom2pipe.manage_composable.data_put') as ad_put_mock, \
                patch('caom2pipe.manage_composable.get_artifact_metadata') as \
                art_mock, \
                patch('caom2pipe.manage_composable.exec_cmd') as exec_mock:
            cadc_client_mock.return_value.data_get.return_value = mc.CadcException(
                'test')
            result = preview_augmentation.visit(obs, **kwargs)
            assert not http_mock.called, 'http mock should not be called'
            assert not ad_put_mock.called, 'ad put mock should not be called'
            assert not art_mock.called, 'art mock should not be called'
            assert not exec_mock.called, 'exec mock should not be called'
            assert result is not None, 'expect a result'
            assert result['artifacts'] == 0, 'no artifacts should be updated'
            assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, \
                'no new artifacts'

        test_rejected.persist_state()
        assert os.path.exists(REJECTED_FILE)
    finally:
        if os.path.exists(REJECTED_FILE):
            os.unlink(REJECTED_FILE)


def test_preview_augment_unknown_no_preview():
    obs = mc.read_obs_from_file(TEST_OBS_FILE)
    obs.planes[TEST_PRODUCT_ID].data_release = datetime.utcnow()
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'

    if os.path.exists(REJECTED_FILE):
        os.unlink(REJECTED_FILE)
    test_rejected = mc.Rejected(REJECTED_FILE)

    cadc_client_mock = Mock()
    kwargs = {'working_directory': TEST_DATA_DIR,
              'cadc_client': cadc_client_mock,
              'stream': 'stream',
              'rejected': test_rejected}

    with patch('caom2pipe.manage_composable.http_get',
               side_effect=mc.CadcException(
                   '404 Client Error: Not Found for url')) as http_mock, \
            patch('caom2pipe.manage_composable.data_put') as ad_put_mock, \
            patch('caom2pipe.manage_composable.get_artifact_metadata') as \
                art_mock, \
            patch('caom2pipe.manage_composable.exec_cmd') as exec_mock:
        cadc_client_mock.return_value.data_get.return_value = mc.CadcException(
            'test')
        try:
            ignore = preview_augmentation.visit(obs, **kwargs)
        except Exception as e:
            pass  # should have got here
        test_url = '{}{}.fits'.format(preview_augmentation.PREVIEW_URL,
                                      TEST_PRODUCT_ID)
        test_prev = '{}/{}.jpg'.format(TEST_DATA_DIR, TEST_PRODUCT_ID)
        http_mock.assert_called_with(test_url, test_prev),  'mock not called'
        assert not ad_put_mock.called, 'ad put mock should not be called'
        assert not art_mock.called, 'art mock should not be called'
        assert not exec_mock.called, 'exec mock should not be called'


# @pytest.mark.skip(reason='Possibly not needed')
# @patch('gem2caom2.GemName._get_obs_id')
# def test_bounds_augment_plane(mock_obs_id):
#     mock_obs_id.return_value = TEST_OBS
#     test_fqn = os.path.join(TESTDATA_DIR, '{}.in.xml'.format(TEST_FP_OBS))
#     test_obs = mc.read_obs_from_file(test_fqn)
#     assert len(test_obs.planes[TEST_FP_OBS].artifacts) == 2
#
#     test_kwargs = {'working_directory': TESTDATA_DIR,
#                    'science_file': TEST_FP_FILE,
#                    'cadc_client': None}
#     test_result = plane_augmentation.visit(test_obs, **test_kwargs)
#
#     # additional Plane metadata for AladinLite
#     # test_obs.planes[TEST_FP_OBS].position.dimension = Dimension2D(3107, 2302)
#     # test_obs.planes[TEST_FP_OBS].position.sample_size = 0.1455995277107371
#
#     validate(test_obs, True)
#
#     # mc.write_obs_to_file(test_obs, '{}-fp.xml'.format(TEST_FP_OBS))
#
#     assert test_result is not None, 'expected update visit return value'
#     assert test_result['planes'] == 1
#
#     assert len(test_obs.planes) == 1
#     assert test_obs.planes[TEST_FP_OBS].position is not None, \
#         'expected Plane.position'
#     assert test_obs.planes[TEST_FP_OBS].position.bounds is not None, \
#         'expected.Plane.position.bounds'
#
#     bounds = test_obs.planes[TEST_FP_OBS].position.bounds
#     assert bounds.samples is not None, 'expected bounds.samples'
#     assert len(bounds.samples.vertices) == 24, 'expected 24 bounds samples'
#     assert bounds.points is not None, 'expected bounds.points'
#     assert len(bounds.points) == 8, 'expected 8 bounds points'
#
#     assert len(test_obs.planes[TEST_FP_OBS].artifacts) == 2
#     fits_uri = 'gemini:GEM/{}'.format(TEST_FP_FILE)
#     artifact = test_obs.planes[TEST_FP_OBS].artifacts[fits_uri]
#     assert artifact.content_type == 'application/fits', \
#         'expected ContentType application/fits'
#     assert artifact.content_length == 15102720, \
#         'expected ContentLength 15102720'
#     assert artifact.content_checksum.checksum == \
#            '5c4f933788ae79553951d10cb5cbedd6', \
#         'expected ContentChecksum 5c4f933788ae79553951d10cb5cbedd6'