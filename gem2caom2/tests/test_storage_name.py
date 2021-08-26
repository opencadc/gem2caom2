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

from os import path
from mock import patch
import gem_mocks
from caom2pipe import manage_composable as mc
from gem2caom2 import GemName, SCHEME, COLLECTION
from gem2caom2 import external_metadata as em


def test_is_valid():
    mock_obs_id = 'GN-2013B-Q-28-150-002'
    assert GemName(file_name='anything.fits', obs_id=mock_obs_id).is_valid()
    assert GemName(file_name='anything.jpg', obs_id=mock_obs_id).is_valid()


@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('caom2pipe.client_composable.query_tap_client')
def test_storage_name(tap_mock, cap_mock):
    tap_mock.side_effect = gem_mocks.mock_query_tap
    cap_mock.return_value = 'https://localhost'
    test_config = mc.Config()
    test_config.proxy_fqn = path.join(gem_mocks.TEST_DATA_DIR, 'cadcproxy.pem')
    test_config.tap_id = 'ivo://cadc.nrc.ca/test'
    em.init_global(test_config)
    mock_obs_id = 'GN-2013B-Q-28-150-002'
    test_sn = GemName(file_name='N20131203S0006i.fits.bz2', obs_id=mock_obs_id)
    assert test_sn.file_uri == f'{SCHEME}:{COLLECTION}/N20131203S0006i.fits'
    assert test_sn.file_name == 'N20131203S0006i.fits'
    assert test_sn.prev == 'N20131203S0006i.jpg'
    assert test_sn.thumb == 'N20131203S0006i_th.jpg'
    assert test_sn.compressed_file_name is None
    assert test_sn.file_id == 'N20131203S0006i'

    test_sn = GemName(file_name='S20060920S0137.jpg', obs_id=mock_obs_id)
    assert test_sn.file_uri == f'{SCHEME}:{COLLECTION}/S20060920S0137.jpg'
    assert test_sn.file_name == 'S20060920S0137.jpg'
    assert test_sn.prev == 'S20060920S0137.jpg'
    assert test_sn.thumb == 'S20060920S0137_th.jpg'
    assert test_sn.compressed_file_name is None

    test_sn = GemName(file_name='N20100104S0208.fits.header')
    assert test_sn.obs_id == 'GN-2009B-Q-121-15-001', 'wrong obs id'
    assert test_sn.file_uri == f'{SCHEME}:{COLLECTION}/N20100104S0208.fits'
    assert (
        test_sn.external_urls
        == 'https://archive.gemini.edu/fullheader/N20100104S0208.fits'
    )

    test_sn = GemName(file_name='N20200810A0490r.fits')
    assert test_sn.obs_id == 'N20200810A0490', 'wrong obs id'
    assert test_sn.product_id == 'N20200810A0490r', 'wrong product id'
    assert test_sn.file_uri == f'{SCHEME}:{COLLECTION}/N20200810A0490r.fits'
    assert (
        test_sn.external_urls
        == 'https://archive.gemini.edu/fullheader/N20200810A0490r.fits'
    )
    assert (
        test_sn.lineage
        == f'{test_sn.obs_id}r/{SCHEME}:{COLLECTION}/{test_sn.file_id}.fits'
    ), 'wrong lineage'

    test_sn = GemName(file_name='SDCH_20200131_0010.fits')
    assert test_sn.obs_id == 'GS-CAL20200131-10-0131', 'wrong obs id'
    assert test_sn.product_id == 'SDC_20200131_0010', 'wrong product id'
    assert test_sn.file_uri == f'{SCHEME}:{COLLECTION}/SDCH_20200131_0010.fits'
    assert (
            test_sn.external_urls
            == 'https://archive.gemini.edu/fullheader/SDCH_20200131_0010.fits'
    )
    assert (
        test_sn.lineage
        == f'SDC_20200131_0010/{SCHEME}:{COLLECTION}/SDCH_20200131_0010.fits'
    ), 'wrong lineage'
