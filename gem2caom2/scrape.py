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

import logging

from astropy.table import Table
from collections import OrderedDict
from datetime import datetime

from caom2pipe import manage_composable as mc


__all__ = ['read_json_file_list_page', 'parse_json_file_list', 'find_direct']

JSON_FILE_LIST = \
    'https://archive.gemini.edu/jsonfilelist/notengineering/NotFail/filepre='
JSON_METADATA = \
    'https://archive.gemini.edu/jsonsummary/canonical/notengineering/NotFail/'


def read_json_file_list_page(start_time_s, last_processed_time_s):
    """

    :param start_time_s: Represents a timestamp for the day being queried.
    :param last_processed_time_s: The archive.gemini.edu timestamp for the
        last file to be processed at CADC. Anything newer should be
        processed.
    :return:
    """
    file_names = {}
    date_str = datetime.fromtimestamp(start_time_s).strftime('%Y%m%d')
    for telescope_str in ['S', 'N']:
        file_list_url = f'{JSON_FILE_LIST}{telescope_str}{date_str}/'
        logging.debug('Querying {}'.format(file_list_url))
        response = None
        try:
            response = mc.query_endpoint(file_list_url)
            if response is None:
                logging.warning(
                    'Could not query {}'.format(file_list_url))
            else:
                temp = parse_json_file_list(
                    response.text, last_processed_time_s)
                temp_file_names = file_names
                file_names = {**temp, **temp_file_names}
                response.close()
        finally:
            if response is not None:
                response.close()
    # order the list of work to be done by time
    ordered_file_names = OrderedDict(file_names.items())
    return ordered_file_names


def parse_json_file_list(json_string, last_processed_time_s):
    work_list = {}
    if not json_string.startswith('[]'):
        temp = Table.read(json_string, format='pandas.json')
        for entry in temp:
            # e.g. 2019-11-01 00:01:34.610517+00:00, and yes, I know about %z
            entry_ts_s = datetime.strptime(
                entry['lastmod'].replace('+00:00', ''),
                '%Y-%m-%d %H:%M:%S.%f').timestamp()
            f_name = entry['filename']
            # the same file name can be in the list returned more than once
            if entry_ts_s >= last_processed_time_s:
                if f_name in work_list.values():
                    for key, value in work_list.items():
                        if f_name == value:
                            existing_ts = key
                            break
                    if entry_ts_s > existing_ts:
                        del work_list[existing_ts]
                        work_list[entry_ts_s] = f_name
                        logging.debug(f'Replacing {entry["filename"]} in work '
                                      f'list.')
                else:
                    work_list[entry_ts_s] = f_name
                    logging.debug(f'Adding {entry["filename"]} to work list.')
    return work_list


def find_direct(data_label):
    metadata = None
    data_label_url = f'{JSON_METADATA}{data_label}'
    logging.debug(f'Querying {data_label_url}')
    response = None
    try:
        response = mc.query_endpoint(data_label_url)
        if response is None:
            logging.warning(f'Could not query {data_label_url}')
        else:
            metadata = response.json()
            response.close()
            if metadata is None or len(metadata) == 0:
                file_url = f'{JSON_METADATA}/filepre={data_label}'
                logging.warning(f'Try by file name {file_url}')
                response = mc.query_endpoint(file_url)
                if response is None:
                    logging.warning(f'Could not query {data_label_url}')
                else:
                    metadata = response.json()
                    response.close()
    finally:
        if response is not None:
            response.close()
    return metadata


def find_data_label_by_file_name(f_name):
    metadata = None
    f_name_url = f'{JSON_METADATA}/filepre={f_name}'
    logging.debug(f'Querying {f_name_url}')
    response = None
    try:
        response = mc.query_endpoint(f_name_url)
        if response is None:
            logging.warning(f'Could not query {f_name_url}')
        else:
            metadata = response.json()
            response.close()
    finally:
        if response is not None:
            response.close()
    return metadata


def parse_for_data_label(json_string, f_name):
    obs_id = None
    last_mod_s = None
    for entry in json_string:
        if entry.get('name') == f_name:
            obs_id = entry.get('data_label')
            last_mod_s = mc.make_seconds(entry.get('lastmod'))
    return obs_id, last_mod_s
