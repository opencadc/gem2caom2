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

from caom2pipe import client_composable as clc
from caom2pipe import data_source_composable as dsc
from caom2pipe import manage_composable as mc
from caom2pipe import reader_composable as rdc
from caom2pipe import run_composable as rc
from gem2caom2 import main_app, preview_augmentation, external_metadata
from gem2caom2 import pull_augmentation, data_source, builder
from gem2caom2 import cleanup_augmentation, fits2caom2_augmentation
from gem2caom2 import gemini_reader

DATA_VISITORS = []
META_VISITORS = [
    fits2caom2_augmentation,
    preview_augmentation,
    pull_augmentation,
    cleanup_augmentation,
]


def _run():
    """
    Uses a todo file with file names, even though Gemini provides
    information about existing data referenced by observation ID.
    """
    config = mc.Config()
    config.get_executors()
    # external_metadata.init_global(config=config)
    clients = None
    if config.use_local_files or mc.TaskType.SCRAPE in config.task_types:
        source = dsc.ListDirSeparateDataSource(config)
        meta_visitors = [
            fits2caom2_augmentation,
            preview_augmentation,
            cleanup_augmentation,
        ]
        metadata_reader = rdc.FileMetadataReader()
    else:
        source = dsc.TodoFileDataSource(config)
        meta_visitors = META_VISITORS
        clients = clc.ClientCollection(config)
        metadata_reader = rdc.StorageClientReader(clients.data_client)
    reader_lookup = gemini_reader.GeminiMetadataLookup(metadata_reader)
    name_builder = builder.GemObsIDBuilder(
        config, metadata_reader, reader_lookup
    )
    return rc.run_by_todo(
        config,
        name_builder,
        chooser=None,
        meta_visitors=meta_visitors,
        source=source,
        metadata_reader=metadata_reader,
        clients=clients,
    )


def run():
    """Wraps _run in exception handling, with sys.exit calls."""
    try:
        result = _run()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_single():
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
        storage_name = builder.GemObsIDBuilder(config).build(sys.argv[1])
    else:
        raise mc.CadcException('No code to handle running GEM by obs id.')
    external_metadata.init_global(config=config)
    return rc.run_single(
        config,
        storage_name,
        main_app.APPLICATION,
        META_VISITORS,
        DATA_VISITORS,
    )


def run_single():
    """
    Run the processing for a single entry.
    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    try:
        result = _run_single()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_by_public():
    """Run the processing for observations that are public, but there are
    no artifacts representing the previews in CAOM, or a FITS file in ad.

    Called as gem_run_public. The time-boxing is based on timestamps from a
    state.yml file. Call once/day, since data release timestamps have times
    of 00:00:00.000.

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    config = mc.Config()
    config.get_executors()
    external_metadata.init_global(config=config)
    name_builder = builder.GemObsIDBuilder(config)
    incremental_source = data_source.PublicIncremental(config)
    return rc.run_by_state(
        config=config,
        name_builder=name_builder,
        bookmark_name=data_source.GEM_BOOKMARK,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        end_time=None,
        source=incremental_source,
        chooser=None,
    )


def run_by_public():
    try:
        result = _run_by_public()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_state():
    """Run incremental processing for observations that are posted on the site
    archive.gemini.edu. TODO in the future this will depend on the incremental
    query endpoint.

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    config = mc.Config()
    config.get_executors()
    state = mc.State(config.state_fqn)
    end_timestamp_s = state.bookmarks.get(data_source.GEM_BOOKMARK).get(
        'end_timestamp', datetime.now()
    )
    end_timestamp_dt = mc.make_time_tz(end_timestamp_s)
    logging.info(f'{main_app.APPLICATION} will end at {end_timestamp_s}')
    # external_metadata.init_global(config=config)
    session = mc.get_endpoint_session()
    metadata_reader = gemini_reader.GeminiMetadataReader(session)
    reader_lookup = gemini_reader.GeminiMetadataLookup(metadata_reader)
    name_builder = builder.GemObsIDBuilder(
        config, metadata_reader, reader_lookup
    )
    incremental_source = data_source.IncrementalSource(
        session, metadata_reader
    )
    result = rc.run_by_state(
        config=config,
        name_builder=name_builder,
        bookmark_name=data_source.GEM_BOOKMARK,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        end_time=end_timestamp_dt,
        source=incremental_source,
        chooser=None,
        metadata_reader=metadata_reader,
    )
    if incremental_source.max_records_encountered:
        logging.warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        logging.warning('Encountered maximum records!!')
        logging.warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        result |= -1
    return result


def run_state():
    try:
        result = _run_state()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
