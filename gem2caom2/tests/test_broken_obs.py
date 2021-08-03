# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
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

import sys

from caom2pipe import manage_composable as mc
from gem2caom2 import main_app, GemName, external_metadata
from gem2caom2.util import Inst, COLLECTION

from mock import patch, Mock
import gem_mocks


@patch('caom2utils.cadc_client_wrapper.StorageClientWrapper')
@patch('gem2caom2.external_metadata.get_obs_metadata')
@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('gem2caom2.program_metadata.get_pi_metadata')
@patch('gem2caom2.external_metadata.DefiningMetadataFinder')
def test_missing_provenance(
    get_mock,
    pi_mock,
    cap_mock,
    get_obs_mock,
    client_mock,
):
    test_config = mc.Config()
    test_config.get_executors()

    external_metadata.init_global(test_config)

    cap_mock.return_value = 'https://localhost'
    pi_mock.side_effect = gem_mocks.mock_get_pi_metadata
    get_mock.return_value.get.side_effect = gem_mocks.mock_get_dm
    get_obs_mock.side_effect = gem_mocks.mock_get_obs_metadata
    client_mock.return_value.info.side_effect = gem_mocks.mock_get_file_info

    test_f_name = 'gS20171114S0185_bias.fits.header'
    test_obs_id = 'GS-CAL20171114-2-086-G-BIAS'
    test_storage_name = GemName(obs_id=test_obs_id, file_name=test_f_name)
    test_fqn = f'{gem_mocks.TEST_DATA_DIR}/broken_files/{test_f_name}'
    actual_fqn = (
        f'{gem_mocks.TEST_DATA_DIR}/broken_files/{test_obs_id}.actual.xml'
    )
    sys.argv = (
        f'{main_app.APPLICATION} --quiet --no_validate '
        f'--local {test_fqn} '
        f'--plugin {gem_mocks.PLUGIN} --module {gem_mocks.PLUGIN} '
        f'--observation {COLLECTION} {test_obs_id} --out {actual_fqn} '
        f'--lineage {test_storage_name.lineage} '
        f'--resource-id ivo://cadc.nrc.ca/test '
    ).split()
    main_app.to_caom2()
    expected_fqn = (
        f'{gem_mocks.TEST_DATA_DIR}/broken_files/{test_obs_id}.expected.xml'
    )
    compare_result = mc.compare_observations(actual_fqn, expected_fqn)
    if compare_result is not None:
        assert False, compare_result
