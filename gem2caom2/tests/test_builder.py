# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
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
import pytest

from mock import Mock, patch

from caom2pipe import manage_composable as mc
from caom2pipe import reader_composable as rdc
from gem2caom2 import SCHEME, V_SCHEME, COLLECTION
from gem2caom2 import builder, gemini_reader
from gem2caom2 import external_metadata as em

import gem_mocks


@patch('gem2caom2.gemini_reader.GeminiMetadataReader._retrieve_headers')
@patch('gem2caom2.gemini_reader.GeminiMetadataReader._retrieve_json')
def test_builder(file_info_mock, header_mock):
    file_info_mock.side_effect = gem_mocks.mock_get_obs_metadata

    test_config = mc.Config()
    test_config.working_directory = '/test_files'
    test_config.proxy_fqn = os.path.join(
        gem_mocks.TEST_DATA_DIR, 'test_proxy.pem'
    )
    session_mock = Mock()
    test_reader = gemini_reader.GeminiMetadataReader(session_mock)
    test_metadata = gemini_reader.GeminiMetadataLookup(test_reader)
    test_subject = builder.GemObsIDBuilder(
        test_config, test_reader, test_metadata
    )

    test_entries = ['S20050825S0143.fits', 'TX20131117_raw.3002.fits']
    for test_entry in test_entries:
        for task_type in [mc.TaskType.INGEST, mc.TaskType.SCRAPE]:
            test_config.task_types = [task_type]
            test_result = test_subject.build(test_entry)
            assert (
                test_result is not None
            ), f'expect a result'
            assert (
                test_result.file_uri
                == f'{SCHEME}:{COLLECTION}/{test_entry}'
            ), 'wrong file uri'
            assert (
                test_result.prev_uri
                == f'{SCHEME}:{COLLECTION}/{test_result.prev}'
            ), 'wrong preview uri'
            assert (
                test_result.thumb_uri
                == f'{V_SCHEME}:{COLLECTION}/{test_result.thumb}'
            ), 'wrong thumb uri'
            assert test_result.obs_id is not None, f'expect an obs id'

    test_config.task_types = [mc.TaskType.INGEST]
    test_entry = 'GN-DATA-LABEL'
    with pytest.raises(mc.CadcException):
        ignore = test_subject.build(test_entry)


@patch('caom2pipe.reader_composable.FileMetadataReader._retrieve_headers')
@patch('caom2pipe.reader_composable.FileMetadataReader._retrieve_file_info')
def test_builder_local(file_info_mock, header_mock):
    file_info_mock.side_effect = gem_mocks.mock_get_obs_metadata
    test_reader = rdc.FileMetadataReader()
    test_metadata = gemini_reader.GeminiMetadataLookup(test_reader)

    try:
        test_config = mc.Config()
        test_config.data_sources = ['/test_files']
        test_config.use_local_files = True
        test_entry = '/test_files/S20191214S0301.fits'
        test_config.task_types = [mc.TaskType.INGEST]
        test_subject = builder.GemObsIDBuilder(
            test_config, test_reader, test_metadata
        )
        test_result = test_subject.build(test_entry)
        assert test_result is not None, 'expect a result'
        assert (
            test_result.file_uri == 'gemini:GEMINI/S20191214S0301.fits'
        ), 'file'
        assert (
            test_result.prev_uri == 'gemini:GEMINI/S20191214S0301.jpg'
        ), 'prev'
        assert (
            test_result.thumb_uri == 'cadc:GEMINI/S20191214S0301_th.jpg'
        ), 'thumb'
        assert (
            test_result.source_names[0] == '/test_files/S20191214S0301.fits'
        ), 'wrong source_names'
    finally:
        em.defining_metadata_finder = None
