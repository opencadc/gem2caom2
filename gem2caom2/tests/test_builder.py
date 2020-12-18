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

from mock import patch

from caom2pipe import manage_composable as mc
from gem2caom2 import A_SCHEME, SCHEME, V_SCHEME, COLLECTION, ARCHIVE
from gem2caom2 import builder
from gem2caom2 import external_metadata as em

import gem_mocks


@patch('gem2caom2.external_metadata.CadcTapClient')
@patch('gem2caom2.external_metadata.get_obs_metadata')
def test_builder(obs_metadata_mock, tap_client_mock):
    obs_metadata_mock.side_effect = gem_mocks.mock_get_obs_metadata

    test_config = mc.Config()
    test_config.working_directory = '/test_files'
    test_config.proxy_fqn = os.path.join(gem_mocks.TEST_DATA_DIR,
                                         'test_proxy.pem')
    em.init_global(incremental=True, config=test_config)
    test_subject = builder.GemObsIDBuilder(test_config)

    test_entry = 'S20050825S0143.fits'
    for support in [False, True]:
        test_config.features.supports_latest_client = support
        test_config.features.use_file_names = True
        for task_type in [mc.TaskType.INGEST, mc.TaskType.SCRAPE]:
            test_config.task_types = [task_type]
            test_result = test_subject.build(test_entry)
            assert test_result is not None, \
                f'expect a result support {support}'
            expected_path = COLLECTION if support else ARCHIVE
            assert test_result.file_uri == \
                   f'{SCHEME}:{expected_path}/{test_entry}', 'wrong file uri'
            assert test_result.prev_uri == \
                   f'{SCHEME}:{expected_path}/{test_result.prev}', \
                   'wrong preview uri'
            expected_scheme = V_SCHEME if support else A_SCHEME
            assert test_result.thumb_uri == \
                   f'{expected_scheme}:{expected_path}/{test_result.thumb}', \
                   'wrong thumb uri'

        test_config.task_types = [mc.TaskType.INGEST]
        test_config.features.use_file_names = False
        with pytest.raises(mc.CadcException):
            test_result = test_subject.build(test_entry)
