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

import os
import pytest

from unittest.mock import patch

import gem_mocks

pytest.main(args=['-s', os.path.abspath(__file__)])

# structured by observation id, list of file ids that make up a multi-plane
# observation
DIR_NAME = 'multi_plane'
LOOKUP = {
    'GS-CAL20101028-5-004': ['mrgS20101028S0134', 'S20101028S0134'],
    'GN-2007B-C-6-5-005': ['TX20071021_RAW.2037', 'TX20071021_SUM.2037'],
    'TX20170321.2505': [
        'TX20170321_raw.2505',
        'TX20170321_red.2505',
        'TX20170321_sum.2505',
    ],
    'TX20170321.2507': ['TX20170321_raw.2507', 'TX20170321_red.2507'],
    'GS-2002B-Q-22-13-0161': [
        '2002dec02_0161',
        'P2002DEC02_0161_SUB.0001',
        'P2002DEC02_0161_SUB',
    ],
    'GN-2015B-Q-1-12-1003': [
        'N20150807G0044',
        'N20150807G0044i',
        'N20150807G0044m',
    ],
    'GN-2012A-Q-124-1-003': ['N20120905S0122', 'N20120905S0122_arc'],
    'GS-2006A-Q-60-11-001': ['S20060412S0056', 'rS20060412S0056'],
    'N20191219A0004': ['N20191219A0004b', 'N20191219A0004r'],
    'N20191219A0001': ['N20191219A0001b', 'N20191219A0001r'],
    'S20200316Z0308': ['S20200316Z0308b', 'S20200316Z0308r'],
    'S20190716Z0925': ['S20190716Z0925b', 'S20190716Z0925r'],
    'N20190313A0002': ['N20190313A0002b', 'N20190313A0002r'],
    'N20200215A0004': ['N20200215A0004b', 'N20200215A0004r'],
    'S20190912Z0264': ['S20190912Z0264b', 'S20190912Z0264r'],
    'GS-2003B-Q-23-17-001': ['S20050916S0159', 'rS20050916S0159'],
    'GS-2013B-Q-75-187-001': ['S20130922S0130', 'S20130922S0130_arc'],
    'GS-2020B-Q-211-52-1119': ['SDCH_20201119_0052', 'SDCK_20201119_0052'],
    'GS-CAL20200202-23-0201': ['SDCH_20200201_0023', 'SDCK_20200201_0023'],
    'GS-CAL20201123-20-1122': ['SDCH_20201122_0020', 'SDCK_20201122_0020'],
}


def pytest_generate_tests(metafunc):
    obs_id_list = []
    for ii in LOOKUP:
        obs_id_list.append(ii)
    metafunc.parametrize('test_name', obs_id_list)


@patch('caom2utils.data_util.get_file_type')
@patch('gem2caom2.gemini_metadata.AbstractGeminiMetadataReader._retrieve_json')
@patch('caom2pipe.astro_composable.get_vo_table_session')
@patch('gem2caom2.gemini_metadata.ProvenanceFinder')
@patch('gem2caom2.program_metadata.get_pi_metadata')
def test_visitor(
    get_pi_mock,
    pf_mock,
    svofps_mock,
    json_mock,
    file_type_mock,
    test_name,
    test_config,
):
    expected_fqn = (
        f'{gem_mocks.TEST_DATA_DIR}/multi_plane/' f'{test_name}.expected.xml'
    )
    test_set = []
    for f_name in LOOKUP[test_name]:
        test_set.append(
            f'{gem_mocks.TEST_DATA_DIR}/multi_plane/{f_name}.fits.header'
        )
    gem_mocks._run_test_common(
        data_sources=[os.path.dirname(test_name)],
        get_pi_mock=get_pi_mock,
        svofps_mock=svofps_mock,
        pf_mock=pf_mock,
        json_mock=json_mock,
        file_type_mock=file_type_mock,
        test_set=test_set,
        expected_fqn=expected_fqn,
        test_config=test_config,
    )
