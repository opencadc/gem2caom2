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

import os
import shutil
from datetime import date
from mock import patch, Mock

from caom2pipe import manage_composable as mc
from gem2caom2 import validator

import gem_mocks


@patch('cadcdata.core.net.BaseWsClient.post')
@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
def test_validator(caps_mock, tap_mock):
    caps_mock.return_value = 'https://sc2.canfar.net/sc2repo'
    tap_response = Mock()
    tap_response.status_code = 200
    tap_response.iter_content.return_value = \
        [b'uri\n'
         b'gemini:GEM/S20170102S0663.fits\n'
         b'gemini:GEM/S20170102S0663.jpg\n'
         b'gemini:GEM/N20120102S0663.fits\n'
         b'gemini:GEM/N20120102S0663.jpg\n'
         b'gemini:GEM/N20191102S0665.fits\n'
         b'gemini:GEM/N20191102S0665.jpg\n']

    ad_response = Mock()
    ad_response.status_code = 200
    ad_response.iter_content.return_value = \
        [b'fileName, ingestDate\n']

    global count
    count = 0

    def _mock_return():
        global count
        if count == 0:
            count = 1
            return tap_response
        else:
            return ad_response

    tap_mock.return_value.__enter__.side_effect = _mock_return

    if not os.path.exists('/usr/src/app/cadcproxy.pem'):
        with open('/usr/src/app/cadcproxy.pem', 'w') as f:
            f.write('proxy content')

    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        test_subject = validator.GeminiValidator()
        test_listing_fqn = \
            f'{test_subject._config.working_directory}/{mc.VALIDATE_OUTPUT}'
        if os.path.exists(test_listing_fqn):
            os.unlink(test_listing_fqn)
        if os.path.exists(test_subject._config.work_fqn):
            os.unlink(test_subject._config.work_fqn)

        test_rejected = f'{gem_mocks.TEST_DATA_DIR}/validate/' \
                        f'test_rejected.yml'
        shutil.copy(test_rejected, test_subject._config.rejected_fqn)

        test_source, test_meta, test_data = test_subject.validate()
        assert test_source is not None, 'expected source result'
        assert len(test_source) == 1037, 'wrong number of source results'
        assert 'rS20111124S0053.fits' in test_source, 'wrong result content'
        assert 'rS20111124S0053.jpg' in test_source, 'wrong result content'

        assert test_meta is not None, 'expected destination result'
        assert len(test_meta) == 2, 'wrong # of destination results'
        assert 'S20170102S0663.fits' in test_meta, \
            'wrong destination content'
        assert os.path.exists(test_listing_fqn), 'should create file record'

        test_subject.write_todo()
        assert os.path.exists(test_subject._config.work_fqn), \
            'should create file record'
        with open(test_subject._config.work_fqn, 'r') as f:
            content = f.readlines()
        content_sorted = sorted(content)
        assert content_sorted[0] == '02jul07.0034.fits\n', 'wrong content'

        # assert False
    finally:
        os.getcwd = getcwd_orig
    assert False


def test_date_file_name():
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        # because randomness in naming
        fnames = {'S20170905S0318.fits': date(2017, 9, 1),
                  'rgS20130103S0098_FRINGE.jpg': date(2013, 1, 1),
                  'GS20141226S0203_BIAS.fits': date(2014, 12, 1),
                  'mrgS20160901S0122_add.jpg': date(2016, 9, 1),
                  'N20160403S0236_flat_pasted.fits': date(2016, 4, 1),
                  'N20141109S0266_BIAS': date(2014, 11, 1),
                  'TX20170321_red.2507.fits': date(2017, 3, 1),
                  'N20170616S0540.fits': date(2017, 6, 1),
                  '02jul07.0186.fits': date(2002, 7, 1),
                  'GN2001BQ013-04.fits': date(2001, 1, 1),
                  '2002APR23_591.fits': date(2002, 4, 1),
                  'r01dec05_007.fits': date(2001, 12, 1),
                  'p2004may20_0048_FLAT.fits': date(2004, 5, 1),
                  'P2003JAN14_0148_DARK.fits': date(2003, 1, 1),
                  'ag2003feb19_6.0001.fits': date(2003, 2, 1),
                  '02jun25.0071.fits': date(2002, 6, 1),
                  '01dec10_1078.fits': date(2001, 12, 1),
                  'c2016may18_sci121.fits': date(2016, 5, 1),
                  '01JUN23_1021.jpg': date(2001, 6, 1),
                  'GS2003BQ031-06.fits': date(2003, 1, 1),
                  'GS2004add003-01.fits': date(2004, 1, 1),
                  'GS2005AQ019-01.fits': date(2005, 1, 1)}
        validate = validator.GeminiValidator()
        for f_name, expected_date in fnames.items():
            result = validate._date_file_name(f_name)
            import logging
            logging.error(f'{result} input {f_name}')
            assert isinstance(result, date), f'{f_name}'
            assert result.year == expected_date.year, f'year fail {f_name}'
            assert result.month == expected_date.month, f'month fail {f_name}'
    finally:
        os.getcwd = getcwd_orig
