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

from gem2caom2.util import Inst
from gem2caom2 import obs_file_relationship

from unittest.mock import patch
import gem_mocks


pytest.main(args=['-s', os.path.abspath(__file__)])


def pytest_generate_tests(metafunc):
    if os.path.exists(gem_mocks.TEST_DATA_DIR):

        file_list = []
        for ii in [
            Inst.GMOS,
            Inst.NIRI,
            Inst.GPI,
            Inst.F2,
            Inst.GSAOI,
            Inst.NICI,
            Inst.TRECS,
            Inst.MICHELLE,
            Inst.GRACES,
            Inst.NIFS,
            Inst.GNIRS,
            Inst.PHOENIX,
            Inst.FLAMINGOS,
            Inst.HRWFS,
            Inst.HOKUPAA,
            Inst.OSCIR,
            Inst.BHROS,
            Inst.CIRPASS,
            Inst.TEXES,
            'processed',
            Inst.ALOPEKE,
            Inst.ZORRO,
            Inst.IGRINS,
        ]:
            walk_dir = _get_inst_name(ii)
            for root, dirs, files in os.walk(
                f'{gem_mocks.TEST_DATA_DIR}/{walk_dir}'
            ):
                for file in files:
                    if file.endswith(".header"):
                        file_list.append(os.path.join(root, file))

        metafunc.parametrize('test_name', file_list)


@patch('caom2utils.data_util.get_file_type')
@patch('gem2caom2.gemini_metadata.AbstractGeminiMetadataReader._retrieve_json')
@patch('caom2pipe.reader_composable.FileMetadataReader._retrieve_headers')
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
    test_name,
):

    test_file_id = obs_file_relationship.remove_extensions(
        os.path.basename(test_name)
    )
    expected_fqn = f'{os.path.dirname(test_name)}/{test_file_id}.expected.xml'

    gem_mocks._run_test_common(
        data_sources=[os.path.dirname(test_name)],
        get_pi_mock=get_pi_mock,
        svofps_mock=svofps_mock,
        headers_mock=headers_mock,
        pf_mock=pf_mock,
        json_mock=json_mock,
        file_type_mock=file_type_mock,
        test_set=[test_name],
        expected_fqn=expected_fqn,
    )


def _get_inst_name(inst):
    walk_dir = inst
    if inst != 'processed' and isinstance(inst, Inst):
        walk_dir = inst.value
    return walk_dir
