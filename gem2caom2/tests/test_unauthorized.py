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

from mock import patch
import gem_mocks


@patch('gem2caom2.gemini_metadata.GeminiFileMetadataReader._retrieve_file_info')
@patch('caom2utils.data_util.get_file_type')
@patch('gem2caom2.gemini_metadata.AbstractGeminiMetadataReader._retrieve_json')
@patch('gem2caom2.gemini_metadata.GeminiFileMetadataReader._retrieve_headers')
@patch('caom2pipe.astro_composable.get_vo_table_session')
@patch('gem2caom2.program_metadata.get_pi_metadata')
@patch('gem2caom2.gemini_metadata.ProvenanceFinder')
def test_visitor(
    pf_mock,
    get_pi_mock,
    svofps_mock,
    headers_mock,
    json_mock,
    file_type_mock,
    file_info_mock,
):
    # test case is unauthorized to retrieve metadata from
    # archive.gemini.edu - so no headers, no file

    file_info_mock.return_value = None
    test_fid = 'S20210518S0022'
    test_f_name = f'{test_fid}.fits'
    # the value of test_fqn doesn't actually matter, there just has to be a
    # non-zero length list
    test_fqn = f'{gem_mocks.TEST_DATA_DIR}/broken_files/{test_f_name}'
    expected_fqn = f'{gem_mocks.TEST_DATA_DIR}/GMOS/{test_fid}.expected.xml'
    gem_mocks._run_test_common(
        data_sources=[f'{gem_mocks.TEST_DATA_DIR}/broken_files'],
        get_pi_mock=get_pi_mock,
        svofps_mock=svofps_mock,
        headers_mock=headers_mock,
        pf_mock=pf_mock,
        json_mock=json_mock,
        file_type_mock=file_type_mock,
        test_set=[test_fqn],
        expected_fqn=expected_fqn,
    )
