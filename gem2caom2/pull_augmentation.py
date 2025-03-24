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

import bz2
import logging
import os
import traceback

from datetime import datetime, timezone
from shutil import copyfileobj

from caom2pipe.astro_composable import check_fitsverify
from caom2pipe.manage_composable import CadcException, http_get, make_datetime
from caom2utils.data_util import get_local_file_headers
from gem2caom2.gemini_metadata import retrieve_headers

FILE_URL = 'https://archive.gemini.edu/file'


class PullVisitor:
    """
    A visitor class that handles the retrieval and storage of files from an external source to CADC storage.

    The archive.gemini.edu harvest is done anonymously, so only retrieve the files if the release date is past.

    This visitor does not update an Observation, although the visit method conforms to the pattern.

    Attributes:
        working_dir (str): The working directory for file operations.
        clients (ClientCollection): The clients required for data operations.
        reporter (ExecutorReporter2): The reporter object for logging and reporting.
        observable (Observable): The observable object from the reporter.
        storage_name (GemName): The storage name object containing file information.
        config (Config): The configuration object.
        result (dict): The result of the pull operation, for observability.
        fqn (str): The fully qualified name of the file on disk for staging between the pull and the put
        json_md5sum (str): The MD5 checksum from the JSON metadata at archive.gemini.edu
        release (datetime): The release date of the file.
        logger (Logger): The logger for the class.

    Methods:
        visit(): Attempts to retrieve and store the file if it is not already at CADC.
        look_pull_and_put(): Checks if a file exists at CADC, retrieves it if not, decompresses, and stores it.
        get_metadata(): Retrieves metadata for the file from the staged copy and removes that copy after
            successful storage. If the file was not retrieved, tries the last resort of /fullheader at
            archive.gemini.edu
    """
    def __init__(self, **kwargs):
        self.working_dir = kwargs.get('working_directory', './')
        self.clients = kwargs.get('clients')
        if self.clients is None:
            logging.warning('Need clients to update. Stopping pull visitor.')
            return
        self.reporter = kwargs.get('reporter')
        if self.reporter is None:
            raise CadcException('Visitor needs a reporter parameter.')
        self.observable = self.reporter._observable
        self.storage_name = kwargs.get('storage_name')
        if self.storage_name is None:
            raise CadcException('Visitor needs a storage_name parameter.')
        self.config = kwargs.get('config')
        if self.config is None:
            raise CadcException('Visitor needs a config parameter.')
        self.result = {'file': 0}
        self.decompressed_fqn = os.path.join(self.working_dir, self.storage_name.file_name).replace('.bz2', '')
        self.json_md5sum = self.storage_name.json_metadata.get(self.storage_name.file_uri).get('data_md5')
        self.release = make_datetime(
            self.storage_name.json_metadata.get(self.storage_name.file_uri).get('release')
        )
        self.logger = logging.getLogger(self.__class__.__name__)

    def visit(self):
        self.logger.debug(f'Begin visit for {self.storage_name.obs_id}.')
        if self.observable.rejected.is_bad_metadata(self.storage_name.obs_id):
            self.logger.info(f'Do not pull the file for {self.storage_name.obs_id} because of bad metadata.')
        elif self.release is None or self.release > datetime.now(tz=timezone.utc).replace(tzinfo=None):
            self.logger.info(f'Do not pull the file for {self.storage_name.obs_id} because it is proprietary.')
            # try /fullheader, in case that much metadata is public
            self.get_metadata()
        else:
            self.look_pull_and_put()
            self.get_metadata()
        self.logger.debug('End visit')
        return None

    def look_pull_and_put(self):
        # this should be a decompressed to decompressed size comparison - the decompressed checksum is available
        # from the /jsonsummary endpoint
        cadc_meta = self.clients.data_client.info(self.storage_name.file_uri)
        if (
            self.json_md5sum is not None and cadc_meta is not None and cadc_meta.md5sum.replace('md5:', '')
            != self.json_md5sum
        ) or cadc_meta is None:
            self.logger.debug(f'Different checksums: Source {self.json_md5sum}, CADC {cadc_meta}')
            url = f'{FILE_URL}/{self.storage_name.file_name}'
            if url.endswith('.fits'):
                url = f'{url}.bz2'
            fqn = os.path.join(self.working_dir, self.storage_name.file_name)
            if fqn.endswith('.fits'):
                fqn = f'{fqn}.bz2'
            try:
                http_get(url, fqn)
                if check_fitsverify(fqn):
                    # PD 10-03-25
                    # decompress before storing, as the bzip2 compression algorithm does not support random access
                    with open(self.decompressed_fqn, 'wb') as f_out, bz2.BZ2File(fqn, 'rb') as f_in:
                        # use shutil to control memory consumption
                        copyfileobj(f_in, f_out)
                    self.clients.data_client.put(os.path.dirname(fqn), self.storage_name.file_uri)
                    self.result = {'file': 1}
                    self.logger.info(f'Retrieved {os.path.basename(fqn)} for storage as {self.storage_name.file_name}')
            except CadcException as e:
                # continue metadata ingestion without the file at CADC
                self.logger.error(e)
                self.logger.debug(traceback.format_exc())
        else:
            self.logger.info(f'{self.storage_name.file_uri} already exists at CADC.')

    def get_metadata(self):
        self.logger.debug('Begin get_metadata')
        if os.path.exists(self.decompressed_fqn):
            if self.storage_name.file_uri not in self.storage_name.metadata:
                self.storage_name.metadata[self.storage_name.file_uri] = get_local_file_headers(
                    self.decompressed_fqn
                )
            self.logger.info(f'Removing {self.decompressed_fqn} after successful storage call.')
            os.unlink(self.decompressed_fqn)

        if self.storage_name.file_uri not in self.storage_name.metadata:
            self.storage_name.metadata[self.storage_name.file_uri] = retrieve_headers(
                self.storage_name, 0, self.logger, self.clients, self.config
            )
        self.logger.debug('End get_metadata')


def visit(observation, **kwargs):
    PullVisitor(**kwargs).visit()
    return observation
