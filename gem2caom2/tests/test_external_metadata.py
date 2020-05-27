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

from astropy.table import Table
from mock import patch

from gem2caom2 import external_metadata as ext_md

import gem_mocks


test_subjects = {
    'H2v=2-1S1_G0220': [ext_md.Inst.NIRI, 'H2S1v2-1-G0220'],
    'Kprime_G0206': [ext_md.Inst.NIRI, 'Kprime-G0206'],
    'H2Oice204_G0242': [ext_md.Inst.NIRI, 'H2Oice2045-G0242'],
    'Jcon(121)_G0232': [ext_md.Inst.NIRI, 'Jcont1207-G0232'],
    'Bra_G0238': [ext_md.Inst.NIRI, 'BrAlpha-G0238'],
    'Bracont_G0237': [ext_md.Inst.NIRI, 'BrAlphaCont-G0237'],
    'Brgamma_G0218': [ext_md.Inst.NIRI, 'BrG-G0218'],
    'Jcon1065_G0239': [ext_md.Inst.NIRI, 'Jcont1065-G0239'],
    'hydrocarb_G0231': [ext_md.Inst.NIRI, 'hydrocarbon-G0231'],
}


def test_repair_filter_name():
    for ii in test_subjects:
        test_result = ext_md._repair_filter_name_for_svo(test_subjects[ii][0],
                                                         ii)
        assert test_result == test_subjects[ii][1], 'wrong value'


@patch('gem2caom2.external_metadata.get_obs_metadata')
@patch('caom2pipe.manage_composable.query_tap_client')
def test_caching_relationship(tap_mock, get_obs_mock):
    ext_md.init_global(incremental=True)
    initial_length = 523
    tap_mock.side_effect = _query_mock_none
    get_obs_mock.side_effect = gem_mocks.mock_get_obs_metadata
    test_subject = ext_md.CachingObsFileRelationship()
    # test an entry that's not in the file, not at CADC, is at
    # archive.gemini.edu
    assert len(test_subject.name_list) == initial_length, 'bad initial length'
    test_result = test_subject.get_obs_id('N20200210S0077.fits')
    assert test_result is not None, 'expect a gemini result'
    assert test_result == 'GN-CAL20200210-22-076', 'wrong gemini result'
    assert len(test_subject.name_list) == initial_length + 1, \
        'bad updated length from Gemini'

    # entry is not in file, but is at CADC
    tap_mock.side_effect = _query_mock_one
    test_result = test_subject.get_obs_id('x.fits')
    assert test_result is not None, 'expect a cadc result'
    assert test_result == 'test_data_label', 'wrong cadc result'
    assert len(test_subject.name_list) == initial_length + 2, \
        'bad updated length from cadc'

    # entry is in file
    test_result = test_subject.get_obs_id('N20170616S0540.fits')
    assert test_result is not None, 'expect a file result'
    assert test_result == 'GN-CAL20170616-11-022', 'wrong file result'
    assert len(test_subject.name_list) == initial_length + 2, \
        'bad updated length from file'


def _query_mock_none(ignore1, ignore2):
    return Table.read('observationID\n'.split('\n'), format='csv')


def _query_mock_one(ignore1, ignore2):
    return Table.read('observationID\n'
                      'test_data_label\n'.split('\n'), format='csv')
