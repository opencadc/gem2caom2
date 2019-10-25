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

import json
import logging
import os

from astropy.io.votable import parse_single_table
from bs4 import BeautifulSoup
from datetime import datetime
from hashlib import md5

from caom2.diff import get_differences
from caom2pipe import manage_composable as mc

import gem2caom2.external_metadata as em
from gem2caom2 import composable


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
FIRST_FILE_LIST = os.path.join(TEST_DATA_DIR, 'S20191010S0.html')
SECOND_FILE_LIST = os.path.join(TEST_DATA_DIR, 'S20191010S3.html')
STATE_FILE = '/usr/src/app/state.yml'
TOO_MANY_FILE_LIST = ''

call_count = 0


def mock_get_votable(url):
    try:
        x = url.split('/')
        filter_name = x[-1].replace('&VERB=0', '')
        votable = parse_single_table(
            f'{TEST_DATA_DIR}/votable/{filter_name}.xml')
        return votable, None
    except Exception as e:
        logging.error(f'get_vo_table failure for url {url}')
        logging.error(e)
        return None, None


def mock_get_pi_metadata(program_id):
    try:
        logging.error(f'program id is {program_id}')
        fname = f'{TEST_DATA_DIR}/programs/{program_id}.xml'
        with open(fname) as f:
            y = f.read()
            soup = BeautifulSoup(y, 'lxml')
            tds = soup.find_all('td')
            if len(tds) > 0:
                title = tds[1].contents[0].replace('\n', ' ')
                pi_name = tds[3].contents[0]
                metadata = {'title': title,
                            'pi_name': pi_name}
                return metadata
        return None
    except Exception as e:
        logging.error(e)
        import traceback
        tb = traceback.format_exc()
        logging.error(tb)


def mock_get_file_info(archive, file_id):
    if '_prev' in file_id:
        return {'size': 10290,
                'md5sum': 'md5:{}'.format(
                    md5('-37'.encode()).hexdigest()),
                'type': 'image/jpeg',
                'name': file_id}
    else:
        return {'size': 665345,
                'md5sum': 'md5:a347f2754ff2fd4b6209e7566637efad',
                'type': 'application/fits',
                'name': file_id}


def mock_get_obs_metadata(file_id):
    try:
        logging.error(f'obs metadata file_id {file_id}')
        fname = f'{TEST_DATA_DIR}/json/{file_id}.json'
        with open(fname) as f:
            y = json.loads(f.read())
            em.om.add(y, file_id)
    except Exception as e:
        logging.error(e)
        import traceback
        tb = traceback.format_exc()
        logging.error(tb)


class Object(object):
    pass

    def close(self):
        pass


def mock_query_endpoint(url, timeout=-1):
    result = Object()
    result.text = None
    global call_count

    file_pre = url.split('/')[-1].replace('filepre=', '')
    if file_pre.startswith('S') and file_pre.endswith('0') and call_count < 2:
        with open(FIRST_FILE_LIST, 'r') as f:
            result.text = f.read()
    elif file_pre.startswith('S') or file_pre.startswith('N'):
        with open(SECOND_FILE_LIST, 'r') as f:
            result.text = f.read()
    elif url == TOO_MANY_FILE_LIST:
        with open(TOO_MANY_FILE_LIST, 'r') as f:
            result.text = f.read()
    else:
        raise Exception('wut {}'.format(url))
    call_count += 1
    return result


def mock_write_state(start_time):
    test_bookmark = {'bookmarks': {composable.GEM_BOOKMARK:
                                   {'last_record': start_time}}}
    mc.write_as_yaml(test_bookmark, STATE_FILE)


def mock_write_state2(prior_timestamp=None):
    # to ensure at least one spin through the execution loop, test case
    # must have a starting time greater than one config.interval prior
    # to 'now', default interval is 10 minutes
    if prior_timestamp is None:
        prior_s = datetime.utcnow().timestamp() - 15 * 60
    else:
        prior_s = mc.make_seconds(prior_timestamp)
    test_start_time = datetime.fromtimestamp(prior_s)
    test_bookmark = {'bookmarks': {composable.GEM_BOOKMARK:
                                   {'last_record': test_start_time}}}
    mc.write_as_yaml(test_bookmark, STATE_FILE)


def mock_repo_create(arg1):
    # arg1 is an Observation instance
    act_fqn = f'{TEST_DATA_DIR}/{arg1.observation_id}.actual.xml'
    ex_fqn = f'{TEST_DATA_DIR}/{arg1.observation_id}.expected.xml'
    mc.write_obs_to_file(arg1, act_fqn)
    result = compare(ex_fqn, act_fqn)
    if result is not None:
        assert False, result


def mock_repo_read(arg1, arg2):
    # arg1 GEMINI arg2 GS-CAL20191010-3-034
    return None


def mock_repo_update():
    return None


def compare(ex_fqn, act_fqn):
    ex = mc.read_obs_from_file(ex_fqn)
    act = mc.read_obs_from_file(act_fqn)
    result = get_differences(ex, act, 'Observation')
    if result:
        result_str = '\n'.join([r for r in result])
        ex_plane = ex.planes.pop()
        msg = f'Differences found obs id {ex.observation_id} ' \
              f'file id {ex_plane.product_id} ' \
              f'instr {ex.instrument.name}\n{result_str}'
        return msg
    return None
