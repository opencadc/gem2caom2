# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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

from unittest.mock import Mock

from gem2caom2 import obs_file_relationship, main_app, util

import gem_mocks


def test_is_processed():
    tests = {
        'c2016may18_sci128': False,
        'abu01aug16_001': False,
        'ag2003feb19_6.0001': True,
        'TX06A_flt.2511': True,
        'GS20141226S0203_BIAS': True,
        'N20070819S0339_dark': True,
        'N20110927S0170_fringe': True,
        'N20120320S0328_stack_fringe': True,
        'N20130404S0512_flat': True,
        'N20140313S0072_flat': True,
        'N20141109S0266_bias': True,
        'N20150804S0348_dark': True,
        'N20160403S0236_flat_pasted': True,
        'S20120922S0406': False,
        'S20131007S0067_fringe': True,
        'S20140124S0039_dark': True,
        'S20141129S0331_dark': True,
        'S20161227S0051': False,
        'fmrgN20020413S0120_add': True,
        'gS20181219S0216_flat': True,
        'gS20190301S0556_bias': True,
        'mfrgS20041117S0073_add': True,
        'mfrgS20160310S0154_add': True,
        'mrgN20041016S0095': True,
        'mrgN20050831S0770_add': True,
        'mrgN20160311S0691_add': True,
        'mrgS20120922S0406': True,
        'mrgS20160901S0122_add': True,
        'mrgS20181016S0184_fringe': True,
        'rS20121030S0136': True,
        'rgS20100212S0301': True,
        'rgS20100316S0366': True,
        'rgS20130103S0098_FRINGE': True,
        'rgS20131109S0166_FRINGE': True,
        'rgS20161227S0051_fringe': True,
        'p2004may20_0048_FLAT': True,
        'p2004may19_0255_COMB': True,
        'P2003JAN14_0148_DARK': True,
        'P2002FEB03_0045_DARK10SEC': True,
        'P2002DEC02_0161_SUB': True,
        'P2002DEC02_0075_SUB.0001': True,
        'N20191219A0004b': False,
        'N20120825S0597_arc': True,
        'rgN20091120S0124_FRINGE': True,
    }
    for ii in tests:
        assert obs_file_relationship.is_processed(ii) == tests[ii], f'failed {ii}'


def test_repair_data_label(test_config):
    # test_config is present so that caom2pipe.StorageName.collection is set properly
    for ii in gem_mocks.LOOKUP.keys():
        test_result = obs_file_relationship.repair_data_label(ii, gem_mocks.LOOKUP[ii][0], gem_mocks.LOOKUP[ii][1])
        if ii == 'S20181230S0025':
            # what happens when an entry is not found
            assert test_result == 'S20181230S0025', (
                f'repair failed for {ii} actual {test_result} expected ' f'{gem_mocks.LOOKUP[ii][0]}'
            )
        elif ii == 'S20201023Z0001b':
            # Alopeke/Zorro have different observationID rules
            assert test_result == 'S20201023Z0001', (
                f'repair failed for {ii} actual {test_result} expected ' f'{gem_mocks.LOOKUP[ii][0]}'
            )
        else:
            assert test_result == gem_mocks.LOOKUP[ii][0], (
                f'repair failed for {ii} actual {test_result} expected ' f'{gem_mocks.LOOKUP[ii][0]}'
            )


def test_repair_data_label_247():
    # unprocessed
    fid1 = 'S20131109S0166'
    dl1 = ['GS-CAL20131109-17-001']

    # processed with a correct data label that's maintained
    fid2 = 'rgS20131109S0166_FRINGE'
    dl2 = ['GS-CAL20131109-17-001-RG-FRINGE']
    #
    # should be GS-CAL20131109-17-001-RG-FRINGE

    # processed with an incorrect data label that has to be fixed
    fid3 = 'rgS20161227S0051_fringe'
    dl3 = ['GS-CAL20161227-5-001', 'GS-CAL20161227-5-001-RG-FRINGE']
    #
    # should be GS-CAL20161227-5-001-RG-FRINGE

    # the unprocessed that goes with the processed
    fid4 = 'S20161227S0051'
    dl4 = ['GS-CAL20161227-5-001']

    # the processed that goes with the unprocessed
    fid5 = 'N20110313S0188_fringe'
    dl5 = ['GN-CAL20110313-900-188']

    d = {
        fid1: dl1,
        fid2: dl2,
        fid3: dl3,
        fid4: dl4,
        fid5: dl5,
    }

    for key, value in d.items():
        index = 0 if len(value) == 1 else 1
        temp = obs_file_relationship.repair_data_label(key, value[0], util.Inst.UNKNOWN)
        assert temp == value[index], f'file id {key} expected {value[index]} actual {temp}'


test_subjects = [
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,1]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,1]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,1]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,1]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,1]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,2]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,2]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,2]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,2]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,2]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,3]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,3]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,3]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,3]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,3]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,4]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,4]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,4]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,4]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,4]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,5]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,5]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,5]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,5]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,5]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,6]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,6]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,6]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,6]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,6]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,7]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,7]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,7]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,7]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,7]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,8]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,8]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,8]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,8]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,8]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,9]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,9]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,9]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,9]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,9]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,10]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,10]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,10]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,10]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,10]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,11]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,11]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,11]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,11]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,11]'],
    ['S20141226S0203', 'tmpfile22889S20141226S0203.fits[SCI,12]'],
    ['S20141226S0204', 'tmpfile22889S20141226S0204.fits[SCI,12]'],
    ['S20141226S0205', 'tmpfile22889S20141226S0205.fits[SCI,12]'],
    ['S20141226S0206', 'tmpfile22889S20141226S0206.fits[SCI,12]'],
    ['S20141226S0207', 'tmpfile22889S20141226S0207.fits[SCI,12]'],
    ['N20070819S0339', 'N20070819S0339.fits[SCI,1]'],
    ['N20070819S0340', 'N20070819S0340.fits[SCI,1]'],
    ['N20070819S0341', 'N20070819S0341.fits[SCI,1]'],
    ['N20070819S0342', 'N20070819S0342.fits[SCI,1]'],
    ['N20070819S0343', 'N20070819S0343.fits[SCI,1]'],
    ['N20070819S0344', 'N20070819S0344.fits[SCI,1]'],
    ['N20070819S0345', 'N20070819S0345.fits[SCI,1]'],
    ['N20070819S0339', 'N20070819S0339.fits[SCI,1]'],
    ['N20070819S0340', 'N20070819S0340.fits[SCI,1]'],
    ['N20070819S0341', 'N20070819S0341.fits[SCI,1]'],
    ['N20070819S0342', 'N20070819S0342.fits[SCI,1]'],
    ['N20070819S0343', 'N20070819S0343.fits[SCI,1]'],
    ['N20070819S0344', 'N20070819S0344.fits[SCI,1]'],
    ['N20070819S0345', 'N20070819S0345.fits[SCI,1]'],
    ['N20070819S0339', 'N20070819S0339.fits[SCI,1]'],
    ['N20070819S0340', 'N20070819S0340.fits[SCI,1]'],
    ['N20070819S0341', 'N20070819S0341.fits[SCI,1]'],
    ['N20070819S0342', 'N20070819S0342.fits[SCI,1]'],
    ['N20070819S0343', 'N20070819S0343.fits[SCI,1]'],
    ['N20070819S0344', 'N20070819S0344.fits[SCI,1]'],
    ['N20070819S0345', 'N20070819S0345.fits[SCI,1]'],
    ['N20130404S0512', 'tmp29851gemcombineN20130404S0512.fits[SCI,1]'],
    ['N20130404S0513', 'tmp29851gemcombineN20130404S0513.fits[SCI,1]'],
    ['N20130404S0514', 'tmp29851gemcombineN20130404S0514.fits[SCI,1]'],
    ['N20130404S0515', 'tmp29851gemcombineN20130404S0515.fits[SCI,1]'],
    ['N20130404S0516', 'tmp29851gemcombineN20130404S0516.fits[SCI,1]'],
    ['N20130404S0517', 'tmp29851gemcombineN20130404S0517.fits[SCI,1]'],
    ['N20130404S0518', 'tmp29851gemcombineN20130404S0518.fits[SCI,1]'],
    ['N20130404S0519', 'tmp29851gemcombineN20130404S0519.fits[SCI,1]'],
    ['N20130404S0520', 'tmp29851gemcombineN20130404S0520.fits[SCI,1]'],
    ['N20130404S0521', 'tmp29851gemcombineN20130404S0521.fits[SCI,1]'],
    ['N20130404S0512', 'tmp29851gemcombineN20130404S0512.fits[SCI,1]'],
    ['N20130404S0513', 'tmp29851gemcombineN20130404S0513.fits[SCI,1]'],
    ['N20130404S0514', 'tmp29851gemcombineN20130404S0514.fits[SCI,1]'],
    ['N20130404S0515', 'tmp29851gemcombineN20130404S0515.fits[SCI,1]'],
    ['N20130404S0516', 'tmp29851gemcombineN20130404S0516.fits[SCI,1]'],
    ['N20130404S0517', 'tmp29851gemcombineN20130404S0517.fits[SCI,1]'],
    ['N20130404S0518', 'tmp29851gemcombineN20130404S0518.fits[SCI,1]'],
    ['N20130404S0519', 'tmp29851gemcombineN20130404S0519.fits[SCI,1]'],
    ['N20130404S0520', 'tmp29851gemcombineN20130404S0520.fits[SCI,1]'],
    ['N20130404S0521', 'tmp29851gemcombineN20130404S0521.fits[SCI,1]'],
    ['N20130404S0512', 'tmp29851gemcombineN20130404S0512.fits[SCI,1]'],
    ['N20130404S0513', 'tmp29851gemcombineN20130404S0513.fits[SCI,1]'],
    ['N20130404S0514', 'tmp29851gemcombineN20130404S0514.fits[SCI,1]'],
    ['N20130404S0515', 'tmp29851gemcombineN20130404S0515.fits[SCI,1]'],
    ['N20130404S0516', 'tmp29851gemcombineN20130404S0516.fits[SCI,1]'],
    ['N20130404S0517', 'tmp29851gemcombineN20130404S0517.fits[SCI,1]'],
    ['N20130404S0518', 'tmp29851gemcombineN20130404S0518.fits[SCI,1]'],
    ['N20130404S0519', 'tmp29851gemcombineN20130404S0519.fits[SCI,1]'],
    ['N20130404S0520', 'tmp29851gemcombineN20130404S0520.fits[SCI,1]'],
    ['N20130404S0521', 'tmp29851gemcombineN20130404S0521.fits[SCI,1]'],
    ['N20141109S0266', 'tmpfile16849_1610gN20141109S0266.fits[SCI,1]'],
    ['N20141109S0267', 'tmpfile16849_1610gN20141109S0267.fits[SCI,1]'],
    ['N20141109S0269', 'tmpfile16849_1610gN20141109S0269.fits[SCI,1]'],
    ['N20141109S0268', 'tmpfile16849_1610gN20141109S0268.fits[SCI,1]'],
    ['N20141109S0270', 'tmpfile16849_1610gN20141109S0270.fits[SCI,1]'],
    ['N20141112S0002', 'tmpfile16849_1610gN20141112S0002.fits[SCI,1]'],
    ['N20141112S0005', 'tmpfile16849_1610gN20141112S0005.fits[SCI,1]'],
    ['N20141112S0001', 'tmpfile16849_1610gN20141112S0001.fits[SCI,1]'],
    ['N20141112S0003', 'tmpfile16849_1610gN20141112S0003.fits[SCI,1]'],
    ['N20141112S0004', 'tmpfile16849_1610gN20141112S0004.fits[SCI,1]'],
    ['N20141112S0093', 'tmpfile16849_1610gN20141112S0093.fits[SCI,1]'],
    ['N20141112S0091', 'tmpfile16849_1610gN20141112S0091.fits[SCI,1]'],
    ['N20141112S0095', 'tmpfile16849_1610gN20141112S0095.fits[SCI,1]'],
    ['N20141112S0092', 'tmpfile16849_1610gN20141112S0092.fits[SCI,1]'],
    ['N20141112S0094', 'tmpfile16849_1610gN20141112S0094.fits[SCI,1]'],
    ['N20141113S0115', 'tmpfile16849_1610gN20141113S0115.fits[SCI,1]'],
    ['N20141113S0116', 'tmpfile16849_1610gN20141113S0116.fits[SCI,1]'],
    ['N20141113S0118', 'tmpfile16849_1610gN20141113S0118.fits[SCI,1]'],
    ['N20141113S0117', 'tmpfile16849_1610gN20141113S0117.fits[SCI,1]'],
    ['N20141113S0119', 'tmpfile16849_1610gN20141113S0119.fits[SCI,1]'],
    ['N20141109S0266', 'tmpfile16849_1610gN20141109S0266.fits[SCI,2]'],
    ['N20141109S0267', 'tmpfile16849_1610gN20141109S0267.fits[SCI,2]'],
    ['N20141109S0269', 'tmpfile16849_1610gN20141109S0269.fits[SCI,2]'],
    ['N20141109S0268', 'tmpfile16849_1610gN20141109S0268.fits[SCI,2]'],
    ['N20141109S0270', 'tmpfile16849_1610gN20141109S0270.fits[SCI,2]'],
    ['N20141112S0002', 'tmpfile16849_1610gN20141112S0002.fits[SCI,2]'],
    ['N20141112S0005', 'tmpfile16849_1610gN20141112S0005.fits[SCI,2]'],
    ['N20141112S0001', 'tmpfile16849_1610gN20141112S0001.fits[SCI,2]'],
    ['N20141112S0003', 'tmpfile16849_1610gN20141112S0003.fits[SCI,2]'],
    ['N20141112S0004', 'tmpfile16849_1610gN20141112S0004.fits[SCI,2]'],
    ['N20141112S0093', 'tmpfile16849_1610gN20141112S0093.fits[SCI,2]'],
    ['N20141112S0091', 'tmpfile16849_1610gN20141112S0091.fits[SCI,2]'],
    ['N20141112S0095', 'tmpfile16849_1610gN20141112S0095.fits[SCI,2]'],
    ['N20141112S0092', 'tmpfile16849_1610gN20141112S0092.fits[SCI,2]'],
    ['N20141112S0094', 'tmpfile16849_1610gN20141112S0094.fits[SCI,2]'],
    ['N20141113S0115', 'tmpfile16849_1610gN20141113S0115.fits[SCI,2]'],
    ['N20141113S0116', 'tmpfile16849_1610gN20141113S0116.fits[SCI,2]'],
    ['N20141113S0118', 'tmpfile16849_1610gN20141113S0118.fits[SCI,2]'],
    ['N20141113S0117', 'tmpfile16849_1610gN20141113S0117.fits[SCI,2]'],
    ['N20141113S0119', 'tmpfile16849_1610gN20141113S0119.fits[SCI,2]'],
    ['N20141109S0266', 'tmpfile16849_1610gN20141109S0266.fits[SCI,3]'],
    ['N20141109S0267', 'tmpfile16849_1610gN20141109S0267.fits[SCI,3]'],
    ['N20141109S0269', 'tmpfile16849_1610gN20141109S0269.fits[SCI,3]'],
    ['N20141109S0268', 'tmpfile16849_1610gN20141109S0268.fits[SCI,3]'],
    ['N20141109S0270', 'tmpfile16849_1610gN20141109S0270.fits[SCI,3]'],
    ['N20141112S0002', 'tmpfile16849_1610gN20141112S0002.fits[SCI,3]'],
    ['N20141112S0005', 'tmpfile16849_1610gN20141112S0005.fits[SCI,3]'],
    ['N20141112S0001', 'tmpfile16849_1610gN20141112S0001.fits[SCI,3]'],
    ['N20141112S0003', 'tmpfile16849_1610gN20141112S0003.fits[SCI,3]'],
    ['N20141112S0004', 'tmpfile16849_1610gN20141112S0004.fits[SCI,3]'],
    ['N20141112S0093', 'tmpfile16849_1610gN20141112S0093.fits[SCI,3]'],
    ['N20141112S0091', 'tmpfile16849_1610gN20141112S0091.fits[SCI,3]'],
    ['N20141112S0095', 'tmpfile16849_1610gN20141112S0095.fits[SCI,3]'],
    ['N20141112S0092', 'tmpfile16849_1610gN20141112S0092.fits[SCI,3]'],
    ['N20141112S0094', 'tmpfile16849_1610gN20141112S0094.fits[SCI,3]'],
    ['N20141113S0115', 'tmpfile16849_1610gN20141113S0115.fits[SCI,3]'],
    ['N20141113S0116', 'tmpfile16849_1610gN20141113S0116.fits[SCI,3]'],
    ['N20141113S0118', 'tmpfile16849_1610gN20141113S0118.fits[SCI,3]'],
    ['N20141113S0117', 'tmpfile16849_1610gN20141113S0117.fits[SCI,3]'],
    ['N20141113S0119', 'tmpfile16849_1610gN20141113S0119.fits[SCI,3]'],
    ['N20141109S0266', 'tmpfile16849_1610gN20141109S0266.fits[SCI,4]'],
    ['N20141109S0267', 'tmpfile16849_1610gN20141109S0267.fits[SCI,4]'],
    ['N20141109S0269', 'tmpfile16849_1610gN20141109S0269.fits[SCI,4]'],
    ['N20141109S0268', 'tmpfile16849_1610gN20141109S0268.fits[SCI,4]'],
    ['N20141109S0270', 'tmpfile16849_1610gN20141109S0270.fits[SCI,4]'],
    ['N20141112S0002', 'tmpfile16849_1610gN20141112S0002.fits[SCI,4]'],
    ['N20141112S0005', 'tmpfile16849_1610gN20141112S0005.fits[SCI,4]'],
    ['N20141112S0001', 'tmpfile16849_1610gN20141112S0001.fits[SCI,4]'],
    ['N20141112S0003', 'tmpfile16849_1610gN20141112S0003.fits[SCI,4]'],
    ['N20141112S0004', 'tmpfile16849_1610gN20141112S0004.fits[SCI,4]'],
    ['N20141112S0093', 'tmpfile16849_1610gN20141112S0093.fits[SCI,4]'],
    ['N20141112S0091', 'tmpfile16849_1610gN20141112S0091.fits[SCI,4]'],
    ['N20141112S0095', 'tmpfile16849_1610gN20141112S0095.fits[SCI,4]'],
    ['N20141112S0092', 'tmpfile16849_1610gN20141112S0092.fits[SCI,4]'],
    ['N20141112S0094', 'tmpfile16849_1610gN20141112S0094.fits[SCI,4]'],
    ['N20141113S0115', 'tmpfile16849_1610gN20141113S0115.fits[SCI,4]'],
    ['N20141113S0116', 'tmpfile16849_1610gN20141113S0116.fits[SCI,4]'],
    ['N20141113S0118', 'tmpfile16849_1610gN20141113S0118.fits[SCI,4]'],
    ['N20141113S0117', 'tmpfile16849_1610gN20141113S0117.fits[SCI,4]'],
    ['N20141113S0119', 'tmpfile16849_1610gN20141113S0119.fits[SCI,4]'],
    ['N20141109S0266', 'tmpfile16849_1610gN20141109S0266.fits[SCI,5]'],
    ['N20141109S0267', 'tmpfile16849_1610gN20141109S0267.fits[SCI,5]'],
    ['N20141109S0269', 'tmpfile16849_1610gN20141109S0269.fits[SCI,5]'],
    ['N20141109S0268', 'tmpfile16849_1610gN20141109S0268.fits[SCI,5]'],
    ['N20141109S0270', 'tmpfile16849_1610gN20141109S0270.fits[SCI,5]'],
    ['N20141112S0002', 'tmpfile16849_1610gN20141112S0002.fits[SCI,5]'],
    ['N20141112S0005', 'tmpfile16849_1610gN20141112S0005.fits[SCI,5]'],
    ['N20141112S0001', 'tmpfile16849_1610gN20141112S0001.fits[SCI,5]'],
    ['N20141112S0003', 'tmpfile16849_1610gN20141112S0003.fits[SCI,5]'],
    ['N20141112S0004', 'tmpfile16849_1610gN20141112S0004.fits[SCI,5]'],
    ['N20141112S0093', 'tmpfile16849_1610gN20141112S0093.fits[SCI,5]'],
    ['N20141112S0091', 'tmpfile16849_1610gN20141112S0091.fits[SCI,5]'],
    ['N20141112S0095', 'tmpfile16849_1610gN20141112S0095.fits[SCI,5]'],
    ['N20141112S0092', 'tmpfile16849_1610gN20141112S0092.fits[SCI,5]'],
    ['N20141112S0094', 'tmpfile16849_1610gN20141112S0094.fits[SCI,5]'],
    ['N20141113S0115', 'tmpfile16849_1610gN20141113S0115.fits[SCI,5]'],
    ['N20141113S0116', 'tmpfile16849_1610gN20141113S0116.fits[SCI,5]'],
    ['N20141113S0118', 'tmpfile16849_1610gN20141113S0118.fits[SCI,5]'],
    ['N20141113S0117', 'tmpfile16849_1610gN20141113S0117.fits[SCI,5]'],
    ['N20141113S0119', 'tmpfile16849_1610gN20141113S0119.fits[SCI,5]'],
    ['N20141109S0266', 'tmpfile16849_1610gN20141109S0266.fits[SCI,6]'],
    ['N20141109S0267', 'tmpfile16849_1610gN20141109S0267.fits[SCI,6]'],
    ['N20141109S0269', 'tmpfile16849_1610gN20141109S0269.fits[SCI,6]'],
    ['N20141109S0268', 'tmpfile16849_1610gN20141109S0268.fits[SCI,6]'],
    ['N20141109S0270', 'tmpfile16849_1610gN20141109S0270.fits[SCI,6]'],
    ['N20141112S0002', 'tmpfile16849_1610gN20141112S0002.fits[SCI,6]'],
    ['N20141112S0005', 'tmpfile16849_1610gN20141112S0005.fits[SCI,6]'],
    ['N20141112S0001', 'tmpfile16849_1610gN20141112S0001.fits[SCI,6]'],
    ['N20141112S0003', 'tmpfile16849_1610gN20141112S0003.fits[SCI,6]'],
    ['N20141112S0004', 'tmpfile16849_1610gN20141112S0004.fits[SCI,6]'],
    ['N20141112S0093', 'tmpfile16849_1610gN20141112S0093.fits[SCI,6]'],
    ['N20141112S0091', 'tmpfile16849_1610gN20141112S0091.fits[SCI,6]'],
    ['N20141112S0095', 'tmpfile16849_1610gN20141112S0095.fits[SCI,6]'],
    ['N20141112S0092', 'tmpfile16849_1610gN20141112S0092.fits[SCI,6]'],
    ['N20141112S0094', 'tmpfile16849_1610gN20141112S0094.fits[SCI,6]'],
    ['N20141113S0115', 'tmpfile16849_1610gN20141113S0115.fits[SCI,6]'],
    ['N20141113S0116', 'tmpfile16849_1610gN20141113S0116.fits[SCI,6]'],
    ['N20141113S0118', 'tmpfile16849_1610gN20141113S0118.fits[SCI,6]'],
    ['N20141113S0117', 'tmpfile16849_1610gN20141113S0117.fits[SCI,6]'],
    ['N20141113S0119', 'tmpfile16849_1610gN20141113S0119.fits[SCI,6]'],
    ['N20150804S0348', 'tmp62119gemcombineN20150804S0348.fits[SCI,1]'],
    ['N20150804S0349', 'tmp62119gemcombineN20150804S0349.fits[SCI,1]'],
    ['N20150804S0350', 'tmp62119gemcombineN20150804S0350.fits[SCI,1]'],
    ['N20150804S0351', 'tmp62119gemcombineN20150804S0351.fits[SCI,1]'],
    ['N20150804S0352', 'tmp62119gemcombineN20150804S0352.fits[SCI,1]'],
    ['N20150804S0353', 'tmp62119gemcombineN20150804S0353.fits[SCI,1]'],
    ['N20150804S0354', 'tmp62119gemcombineN20150804S0354.fits[SCI,1]'],
    ['N20150804S0355', 'tmp62119gemcombineN20150804S0355.fits[SCI,1]'],
    ['N20150804S0356', 'tmp62119gemcombineN20150804S0356.fits[SCI,1]'],
    ['N20150804S0357', 'tmp62119gemcombineN20150804S0357.fits[SCI,1]'],
    ['N20150804S0348', 'tmp62119gemcombineN20150804S0348.fits[SCI,1]'],
    ['N20150804S0349', 'tmp62119gemcombineN20150804S0349.fits[SCI,1]'],
    ['N20150804S0350', 'tmp62119gemcombineN20150804S0350.fits[SCI,1]'],
    ['N20150804S0351', 'tmp62119gemcombineN20150804S0351.fits[SCI,1]'],
    ['N20150804S0352', 'tmp62119gemcombineN20150804S0352.fits[SCI,1]'],
    ['N20150804S0353', 'tmp62119gemcombineN20150804S0353.fits[SCI,1]'],
    ['N20150804S0354', 'tmp62119gemcombineN20150804S0354.fits[SCI,1]'],
    ['N20150804S0355', 'tmp62119gemcombineN20150804S0355.fits[SCI,1]'],
    ['N20150804S0356', 'tmp62119gemcombineN20150804S0356.fits[SCI,1]'],
    ['N20150804S0357', 'tmp62119gemcombineN20150804S0357.fits[SCI,1]'],
    ['N20150804S0348', 'tmp62119gemcombineN20150804S0348.fits[SCI,1]'],
    ['N20150804S0349', 'tmp62119gemcombineN20150804S0349.fits[SCI,1]'],
    ['N20150804S0350', 'tmp62119gemcombineN20150804S0350.fits[SCI,1]'],
    ['N20150804S0351', 'tmp62119gemcombineN20150804S0351.fits[SCI,1]'],
    ['N20150804S0352', 'tmp62119gemcombineN20150804S0352.fits[SCI,1]'],
    ['N20150804S0353', 'tmp62119gemcombineN20150804S0353.fits[SCI,1]'],
    ['N20150804S0354', 'tmp62119gemcombineN20150804S0354.fits[SCI,1]'],
    ['N20150804S0355', 'tmp62119gemcombineN20150804S0355.fits[SCI,1]'],
    ['N20150804S0356', 'tmp62119gemcombineN20150804S0356.fits[SCI,1]'],
    ['N20150804S0357', 'tmp62119gemcombineN20150804S0357.fits[SCI,1]'],
    ['N20160403S0236', 'rgN20160403S0236.fits[SCI,2]'],
    ['N20160403S0237', 'rgN20160403S0237.fits[SCI,2]'],
    ['N20160403S0238', 'rgN20160403S0238.fits[SCI,2]'],
    ['N20160403S0235', 'rgN20160403S0235.fits[SCI,2]'],
    ['N20160403S0240', 'rgN20160403S0240.fits[SCI,2]'],
    ['N20160403S0239', 'rgN20160403S0239.fits[SCI,2]'],
    ['N20160403S0234', 'rgN20160403S0234.fits[SCI,2]'],
    ['N20160403S0241', 'rgN20160403S0241.fits[SCI,2]'],
    ['N20160403S0229', 'rgN20160403S0229.fits[SCI,2]'],
    ['N20160403S0228', 'rgN20160403S0228.fits[SCI,2]'],
    ['N20160403S0230', 'rgN20160403S0230.fits[SCI,2]'],
    ['N20160403S0231', 'rgN20160403S0231.fits[SCI,2]'],
    ['N20160403S0233', 'rgN20160403S0233.fits[SCI,2]'],
    ['N20160403S0232', 'rgN20160403S0232.fits[SCI,2]'],
    ['N20160404S0141', 'rgN20160404S0141.fits[SCI,2]'],
    ['N20160404S0140', 'rgN20160404S0140.fits[SCI,2]'],
    ['N20160404S0139', 'rgN20160404S0139.fits[SCI,2]'],
    ['N20160404S0142', 'rgN20160404S0142.fits[SCI,2]'],
    ['N20160404S0143', 'rgN20160404S0143.fits[SCI,2]'],
    ['N20160404S0145', 'rgN20160404S0145.fits[SCI,2]'],
    ['N20160404S0144', 'rgN20160404S0144.fits[SCI,2]'],
    ['N20160404S0138', 'rgN20160404S0138.fits[SCI,2]'],
    ['N20160404S0137', 'rgN20160404S0137.fits[SCI,2]'],
    ['N20160404S0136', 'rgN20160404S0136.fits[SCI,2]'],
    ['N20160404S0135', 'rgN20160404S0135.fits[SCI,2]'],
    ['S20131007S0067', 'tmp71808gemcombineS20131007S0067.fits[SCI,1]'],
    ['S20131007S0068', 'tmp71808gemcombineS20131007S0068.fits[SCI,1]'],
    ['S20131007S0069', 'tmp71808gemcombineS20131007S0069.fits[SCI,1]'],
    ['S20131007S0067', 'tmp71808gemcombineS20131007S0067.fits[SCI,1]'],
    ['S20131007S0068', 'tmp71808gemcombineS20131007S0068.fits[SCI,1]'],
    ['S20131007S0069', 'tmp71808gemcombineS20131007S0069.fits[SCI,1]'],
    ['S20131007S0067', 'tmp71808gemcombineS20131007S0067.fits[SCI,1]'],
    ['S20131007S0068', 'tmp71808gemcombineS20131007S0068.fits[SCI,1]'],
    ['S20131007S0069', 'tmp71808gemcombineS20131007S0069.fits[SCI,1]'],
    ['S20140124S0039', 'tmp67553gemcombineS20140124S0039.fits[SCI,1]'],
    ['S20140124S0040', 'tmp67553gemcombineS20140124S0040.fits[SCI,1]'],
    ['S20140124S0041', 'tmp67553gemcombineS20140124S0041.fits[SCI,1]'],
    ['S20140124S0042', 'tmp67553gemcombineS20140124S0042.fits[SCI,1]'],
    ['S20140124S0043', 'tmp67553gemcombineS20140124S0043.fits[SCI,1]'],
    ['S20140124S0044', 'tmp67553gemcombineS20140124S0044.fits[SCI,1]'],
    ['S20140124S0039', 'tmp67553gemcombineS20140124S0039.fits[SCI,1]'],
    ['S20140124S0040', 'tmp67553gemcombineS20140124S0040.fits[SCI,1]'],
    ['S20140124S0041', 'tmp67553gemcombineS20140124S0041.fits[SCI,1]'],
    ['S20140124S0042', 'tmp67553gemcombineS20140124S0042.fits[SCI,1]'],
    ['S20140124S0043', 'tmp67553gemcombineS20140124S0043.fits[SCI,1]'],
    ['S20140124S0044', 'tmp67553gemcombineS20140124S0044.fits[SCI,1]'],
    ['S20140124S0039', 'tmp67553gemcombineS20140124S0039.fits[SCI,1]'],
    ['S20140124S0040', 'tmp67553gemcombineS20140124S0040.fits[SCI,1]'],
    ['S20140124S0041', 'tmp67553gemcombineS20140124S0041.fits[SCI,1]'],
    ['S20140124S0042', 'tmp67553gemcombineS20140124S0042.fits[SCI,1]'],
    ['S20140124S0043', 'tmp67553gemcombineS20140124S0043.fits[SCI,1]'],
    ['S20140124S0044', 'tmp67553gemcombineS20140124S0044.fits[SCI,1]'],
    ['S20141129S0331', 'tmp5862gemcombineS20141129S0331.fits[SCI,1]'],
    ['S20141129S0332', 'tmp5862gemcombineS20141129S0332.fits[SCI,1]'],
    ['S20141129S0333', 'tmp5862gemcombineS20141129S0333.fits[SCI,1]'],
    ['S20141129S0334', 'tmp5862gemcombineS20141129S0334.fits[SCI,1]'],
    ['S20141129S0335', 'tmp5862gemcombineS20141129S0335.fits[SCI,1]'],
    ['S20141129S0336', 'tmp5862gemcombineS20141129S0336.fits[SCI,1]'],
    ['S20141129S0337', 'tmp5862gemcombineS20141129S0337.fits[SCI,1]'],
    ['S20141129S0331', 'tmp5862gemcombineS20141129S0331.fits[SCI,1]'],
    ['S20141129S0332', 'tmp5862gemcombineS20141129S0332.fits[SCI,1]'],
    ['S20141129S0333', 'tmp5862gemcombineS20141129S0333.fits[SCI,1]'],
    ['S20141129S0334', 'tmp5862gemcombineS20141129S0334.fits[SCI,1]'],
    ['S20141129S0335', 'tmp5862gemcombineS20141129S0335.fits[SCI,1]'],
    ['S20141129S0336', 'tmp5862gemcombineS20141129S0336.fits[SCI,1]'],
    ['S20141129S0337', 'tmp5862gemcombineS20141129S0337.fits[SCI,1]'],
    ['S20141129S0331', 'tmp5862gemcombineS20141129S0331.fits[SCI,1]'],
    ['S20141129S0332', 'tmp5862gemcombineS20141129S0332.fits[SCI,1]'],
    ['S20141129S0333', 'tmp5862gemcombineS20141129S0333.fits[SCI,1]'],
    ['S20141129S0334', 'tmp5862gemcombineS20141129S0334.fits[SCI,1]'],
    ['S20141129S0335', 'tmp5862gemcombineS20141129S0335.fits[SCI,1]'],
    ['S20141129S0336', 'tmp5862gemcombineS20141129S0336.fits[SCI,1]'],
    ['S20141129S0337', 'tmp5862gemcombineS20141129S0337.fits[SCI,1]'],
    ['S20181219S0216', 'rgS20181219S0216[SCI,1]'],
    ['S20190301S0556', 'tmpfile31966S20190301S0556.fits[SCI,1]'],
    ['S20041117S0073', 'mfrgS20041117S0073_trn'],
    ['S20041117S0074', 'mfrgS20041117S0074_trn'],
    ['S20160310S0154', 'mfrgS20160310S0154_trn'],
    ['S20160310S0155', 'mfrgS20160310S0155_trn'],
    ['S20160310S0156', 'mfrgS20160310S0156_trn'],
    ['S20160310S0157', 'mfrgS20160310S0157_trn'],
    ['S20160310S0158', 'mfrgS20160310S0158_trn'],
    ['S20160310S0160', 'mfrgS20160310S0160_trn'],
    ['N20160311S0691', 'mrgN20160311S0691_trn'],
    ['N20160311S0692', 'mrgN20160311S0692_trn'],
    ['N20160311S0693', 'mrgN20160311S0693_trn'],
    ['N20160311S0694', 'mrgN20160311S0694_trn'],
    ['S20160901S0122', 'mrgS20160901S0122_trn'],
    ['S20160901S0123', 'mrgS20160901S0123_trn'],
    ['S20160901S0124', 'mrgS20160901S0124_trn'],
    ['S20160901S0125', 'mrgS20160901S0125_trn'],
    ['2004may20_0048', 'rawdir$2004may20_0048.fits'],
    ['2004may20_0049', 'rawdir$2004may20_0049.fits'],
    ['2004may20_0050', 'rawdir$2004may20_0050.fits'],
    ['N20060130S0149', 'mrgN20060130S0149_trn'],
]


def test_repair_provenance():
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=gem_mocks.TEST_DATA_DIR)
    try:
        for ii in test_subjects:
            temp = main_app._remove_processing_detritus([ii[1]], 'log_value')
            assert temp[0] == ii[0], f'error {temp[0]} should be {ii[0]}'
    finally:
        os.getcwd = getcwd_orig


def test_repair_data_label_2():
    repairs = {
        'rgS20180122S0236_fringe.fits': [
            'GS-CAL20180122-1-001',
            'GS-CAL20180122-1-001-RG-FRINGE',
        ],
        'gS20150906S0307_bias.fits': [
            'GS-CAL20150906-1-086-BIAS/MBIAS/G-BIAS',
            'GS-CAL20150906-1-086-g-bias',
        ],
    }
    for f_name in repairs.keys():
        result = obs_file_relationship.repair_data_label(f_name, repairs[f_name][0], util.Inst.UNKNOWN)
        assert repairs[f_name][1] == result, f'{result} should have been {repairs[f_name][1]}'
