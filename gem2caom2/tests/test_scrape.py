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

from datetime import datetime
from mock import patch

from caom2pipe import manage_composable as mc

from gem2caom2 import scrape

import gem_mocks as gm

DATA_FILE = f'{gm.TEST_DATA_DIR}/jsonfilelist.json'


def test_parse_json_file_list():
    with open(DATA_FILE, 'r') as f:
        json_string = f.read()

    test_file_name = 'N20191101S0376.fits'
    test_time_end = _get_end_time()
    result = scrape.parse_json_file_list(json_string, test_time_end)
    assert result is not None, 'expected result'
    assert len(result) == 346, 'wrong number of results'
    assert result[test_time_end] == test_file_name, 'wrong entries'


def test_read_json_file_list_page():
    test_file_name = 'N20191101S0119.fits'
    test_end_time = _get_end_time()
    test_start_time = test_end_time - 100000

    with patch('caom2pipe.manage_composable.query_endpoint') as query_mock:
        query_mock.side_effect = _mock_endpoint
        test_work_list = scrape.read_json_file_list_page(
            test_start_time,
            test_end_time)
        assert test_work_list is not None, 'expected result'
        assert len(test_work_list) == 346
        first_entry = test_work_list.popitem(last=False)
        assert first_entry[1] == test_file_name


def _get_end_time():
    return datetime.strptime('2019-11-01 18:22:48.477300',
                             '%Y-%m-%d %H:%M:%S.%f').timestamp()


def _mock_endpoint(url, timeout=-1):
    result = gm.Object()
    result.text = None

    if url.startswith(scrape.JSON_FILE_LIST):
        with open(DATA_FILE, 'r') as f:
            result.text = f.read()
    else:
        raise mc.CadcException('wut?')
    return result
