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
import traceback

from os import path

from caom2pipe import manage_composable as mc
from caom2pipe import name_builder_composable as nbc
from gem2caom2 import gem_name
from gem2caom2.util import COLLECTION, SCHEME
from gem2caom2 import external_metadata as em


__all__ = ['GemObsIDBuilder']


class GemObsIDBuilder(nbc.StorageNameBuilder):
    """
    To be able to build a StorageName instance with an observation ID.
    """

    def __init__(self, config, metadata_reader, metadata):
        super(GemObsIDBuilder, self).__init__()
        self._config = config
        self._metadata_reader = metadata_reader
        self._metadata = metadata
        self._metadata.reader = self._metadata_reader
        self._logger = logging.getLogger(self.__class__.__name__)

    def build(self, entry):
        """
        :param entry: str a Gemini file name, or a fully-qualified file
            name on disk.
        :return: an instance of StorageName for use in execute_composable.
        """
        self._logger.debug(f'Build a StorageName instance for {entry}.')
        try:
            f_name = entry
            if entry != path.basename(entry):
                f_name = path.basename(entry)
            # uri = mc.build_uri(COLLECTION, f_name, SCHEME)
            # metadata = em.defining_metadata_finder.get(uri)
            if (
                mc.TaskType.SCRAPE in self._config.task_types
                or self._config.use_local_files
            ):
                result = gem_name.GemName(
                    file_name=f_name,
                    entry=entry,
                )
                result.source_names = [entry]
            elif '.fits' in entry or '.jpg' in entry:
                result = gem_name.GemName(
                    file_name=f_name,
                    entry=entry,
                )
                result.source_names = [result.file_id]
            else:
                raise mc.CadcException(
                    'The need has not been encountered in the real world '
                    'yet.'
                )
            self._metadata_reader.set(result)
            result.obs_id = self._metadata.data_label
            # result.obs_id = self._metadata_reader.headers.get('DATALAB')
            # if (
            #     result.obs_id is None
            #     and hasattr(self._metadata_reader, 'json_metadata')
            # ):
            #     result.obs_id = self._metadata_reader.json_metadata.get(
            #         'data_label'
            #     )
            self._logger.debug('Done build.')
            return result
        except Exception as e:
            self._logger.error(e)
            self._logger.debug(traceback.format_exc())
            raise mc.CadcException(e)
