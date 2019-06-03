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
import sys
import tempfile
import traceback

from datetime import datetime

from caom2pipe import execute_composable as ec
from caom2pipe import manage_composable as mc
from gem2caom2 import APPLICATION, work
from gem2caom2 import preview_augmentation
from gem2caom2.gem_name import GemName, COLLECTION, ARCHIVE


meta_visitors = [preview_augmentation]
data_visitors = []


def run():
    proxy = '/usr/src/app/cadcproxy.pem'
    ec.run_by_file(GemName, APPLICATION, COLLECTION, proxy,
                   meta_visitors, data_visitors, archive=ARCHIVE)


def run_proxy():
    proxy = '/usr/src/app/cadcproxy.pem'
    ec.run_by_file(GemName, APPLICATION, COLLECTION, proxy,
                   meta_visitors, data_visitors, archive=ARCHIVE)


def run_single():
    """
    Run the processing for a single entry.
    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    config = mc.Config()
    config.get_executors()
    config.resource_id = 'ivo://cadc.nrc.ca/sc2repo'
    if config.features.run_in_airflow:
        temp = tempfile.NamedTemporaryFile()
        mc.write_to_file(temp.name, sys.argv[2])
        config.proxy = temp.name
    else:
        config.proxy = sys.argv[2]
    config.stream = 'default'
    if config.features.use_file_names:
        storage_name = GemName(file_name=sys.argv[1])
    else:
        raise mc.CadcException('No code to handle running GEM by obs id.')
    result = ec.run_single(config, storage_name, APPLICATION, meta_visitors,
                           data_visitors)
    sys.exit(result)


def run_query_2():
    """
    Run the processing for all the entries returned from a time-boxed ad
    query.

    :param sys.argv[1] the timestamp for the > comparison in the time-boxed
        query
    :param sys.argv[2] the timestamp for the <= comparison in the time-boxed
        query

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    prev_exec_date = sys.argv[1]
    exec_date = sys.argv[2]

    config = mc.Config()
    config.get()
    config.stream = 'default'

    file_list = mc.read_file_list_from_archive(config, APPLICATION,
                                               prev_exec_date, exec_date)
    sys.argv = sys.argv[:1]
    result = 0
    if len(file_list) > 0:
        mc.write_to_file(config.work_fqn, '\n'.join(file_list))
        result |= ec.run_by_file(GemName, APPLICATION, COLLECTION,
                                 config.proxy_fqn, meta_visitors,
                                 data_visitors, archive=ARCHIVE)
    sys.exit(result)


def _run_query():
    """
    Run the processing for all the entries returned from a time-boxed CAOM
    query.

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """

    config = mc.Config()
    config.get()
    state = mc.State(config.state_fqn)
    start_time = state.get_bookmark('gemini_timestamp')

    logger = logging.getLogger()
    logger.setLevel(config.logging_level)

    prev_exec_date = start_time
    exec_date = start_time
    now_s = datetime.utcnow()

    logging.debug('Starting at {}'.format(start_time))

    result = 0
    count = 0
    while exec_date < now_s:
        logging.debug(
            'Processing from {} to {}'.format(prev_exec_date, exec_date))
        obs_ids = work.read_obs_ids_from_caom(config, prev_exec_date, exec_date)
        if len(obs_ids) > 0:
            mc.write_to_file(config.work_fqn, '\n'.join(obs_ids))
            result |= ec.run_by_file(GemName, APPLICATION, COLLECTION,
                                     config.proxy_fqn, meta_visitors,
                                     data_visitors, archive=ARCHIVE)
            count += 1
            if count % config.interval == 0:
                state.save_state('gemini_timestamp', prev_exec_date)
                logging.info('Saving timestamp {}'.format(prev_exec_date))
        prev_exec_date = exec_date
        exec_date = mc.increment_time(prev_exec_date, config.interval)

    state.save_state('gemini_timestamp', prev_exec_date)
    logging.info(
        'Done {}, saved state is {}'.format(APPLICATION, prev_exec_date))

    return result


def run_query():
    try:
        result = _run_query()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
