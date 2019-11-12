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
import logging
import os
import sys

import pytest

from shutil import copyfile

import gem2caom2.external_metadata as em

from gem2caom2 import main_app, gem_name
from caom2pipe import manage_composable as mc

from mock import patch, Mock
import gem_mocks


pytest.main(args=['-s', os.path.abspath(__file__)])


def pytest_generate_tests(metafunc):
    if os.path.exists(gem_mocks.TEST_DATA_DIR):

        file_list = []
        for ii in [em.Inst.GMOS, em.Inst.NIRI, em.Inst.GPI, em.Inst.F2,
                   em.Inst.GSAOI, em.Inst.NICI, em.Inst.TRECS, em.Inst.MICHELLE,
                   em.Inst.GRACES, em.Inst.NIFS, em.Inst.GNIRS, em.Inst.PHOENIX,
                   em.Inst.FLAMINGOS, em.Inst.HRWFS, em.Inst.HOKUPAA,
                   em.Inst.OSCIR, em.Inst.BHROS, em.Inst.CIRPASS, em.Inst.TEXES,
                   'processed']:
            walk_dir = _get_inst_name(ii)
            for root, dirs, files in os.walk(
                    '{}/{}'.format(gem_mocks.TEST_DATA_DIR, walk_dir)):
                for file in files:
                    if file.endswith(".header"):
                        file_list.append(os.path.join(root, file))

        metafunc.parametrize('test_name', file_list)


@patch('sys.exit', Mock())
def test_main_app(test_name):
    em.set_ofr(None)
    test_data_size = os.stat(
        os.path.join(gem_mocks.TEST_DATA_DIR, 'from_paul.txt'))
    app_size = os.stat('/app/data/from_paul.txt')
    if test_data_size.st_size != app_size.st_size:
        copyfile(os.path.join(gem_mocks.TEST_DATA_DIR, 'from_paul.txt'),
                 '/app/data/from_paul.txt')
    basename = os.path.basename(test_name)
    dirname = os.path.dirname(test_name)
    file_id = _get_file_id(basename)
    obs_id = _get_obs_id(file_id)
    product_id = file_id
    lineage = _get_lineage(dirname, basename, product_id, file_id)
    input_file = '{}.in.xml'.format(product_id)
    actual_fqn = _get_actual_file_name(dirname, product_id)

    local = _get_local(test_name)
    plugin = gem_mocks.PLUGIN

    with patch('caom2utils.fits2caom2.CadcDataClient') as data_client_mock, \
        patch('gem2caom2.external_metadata.get_obs_metadata') as \
            gemini_client_mock, \
            patch('gem2caom2.external_metadata.get_pi_metadata') as \
            gemini_pi_mock, \
            patch('gem2caom2.svofps.get_vo_table') as svofps_mock:

        data_client_mock.return_value.get_file_info.side_effect = \
            gem_mocks.mock_get_file_info
        gemini_client_mock.side_effect = gem_mocks.mock_get_obs_metadata
        gemini_pi_mock.side_effect = gem_mocks.mock_get_pi_metadata
        svofps_mock.side_effect = gem_mocks.mock_get_votable

        if os.path.exists(actual_fqn):
            os.remove(actual_fqn)

        if os.path.exists(os.path.join(dirname, input_file)):
            sys.argv = \
                ('{} --quiet --no_validate --local {} '
                 '--plugin {} --module {} --in {}/{} --out {} --lineage {}'.
                 format(main_app.APPLICATION, local, plugin, plugin, dirname,
                        input_file, actual_fqn, lineage)).split()
        else:
            sys.argv = \
                ('{} --quiet --no_validate --local {} '
                 '--plugin {} --module {} --observation {} {} --out {} '
                 '--lineage {}'.
                 format(main_app.APPLICATION, local, plugin, plugin,
                        main_app.COLLECTION, obs_id, actual_fqn,
                        lineage)).split()
        print(sys.argv)
        main_app.main_app2()
        expected_fqn = _get_expected_file_name(dirname, product_id)
        compare_result = gem_mocks.compare(
            expected_fqn, actual_fqn, product_id)
        assert compare_result is None, 'compare fail'
        # assert False  # cause I want to see logging messages


def _get_obs_id(file_id):
    return gem_mocks.LOOKUP[file_id][0]


def _get_instr(file_id):
    temp = gem_mocks.LOOKUP[file_id][1]
    return _get_inst_name(temp)


def _get_program_id(file_id):
    return gem_mocks.LOOKUP[file_id][2]


def _get_local(test_name):
    jpg = test_name.replace('.fits.header', '.jpg')
    header_name = test_name
    if os.path.exists(jpg):
        return '{} {}'.format(jpg, header_name)
    else:
        return header_name


def _get_file_id(basename):
    if basename.endswith('jpg'):
        return basename.split('.jpg')[0]
    else:
        return basename.split('.fits')[0]


def _get_lineage(dirname, basename, product_id, file_id):
    logging.error('basename is {}'.format(basename))
    jpg_file = basename.replace('.fits.header', '.jpg')
    if os.path.exists(os.path.join(dirname, jpg_file)):
        jpg = mc.get_lineage(gem_name.ARCHIVE, product_id,
                             '{}.jpg'.format(file_id), gem_name.SCHEME)
        fits = mc.get_lineage(gem_name.ARCHIVE, product_id,
                              '{}.fits'.format(file_id), gem_name.SCHEME)
        return '{} {}'.format(jpg, fits)
    else:
        return mc.get_lineage(gem_name.ARCHIVE, product_id,
                              '{}.fits'.format(file_id), gem_name.SCHEME)


def _get_expected_file_name(dirname, product_id):
    return '{}/{}.xml'.format(dirname, product_id)


def _get_actual_file_name(dirname, product_id):
    return '{}/{}.actual.xml'.format(dirname, product_id)


def _get_inst_name(inst):
    walk_dir = inst
    if inst != 'processed' and isinstance(inst, em.Inst):
        walk_dir = inst.value
    return walk_dir
