# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2022.                            (c) 2022.
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
#  : 4 $
#
# ***********************************************************************
#

import requests

from cadcdata import FileInfo
from caom2utils import data_util
from caom2pipe import manage_composable as mc
from caom2pipe import reader_composable as rdc
from gem2caom2 import obs_file_relationship


__all__ = [
    'GeminiMetadataLookup',
    'GeminiMetadataReader',
    'GEMINI_METADATA_URL',
    'HEADER_URL',
]


GEMINI_METADATA_URL = (
    'https://archive.gemini.edu/jsonsummary/canonical/filepre='
)
HEADER_URL = 'https://archive.gemini.edu/fullheader/'


class GeminiMetadataReader(rdc.MetadataReader):

    def __init__(self, http_session):
        super().__init__()
        self._json_metadata = {}
        self._session = http_session

    @property
    def json_metadata(self):
        return self._json_metadata

    def _retrieve_file_info(self, source_name):
        pass

    def _retrieve_headers(self, source_name):
        self._logger.debug(f'Begin _retrieve_headers for {source_name}')
        header_url = f'{HEADER_URL}{source_name}.fits'
        # Open the URL and fetch the JSON document for the observation
        response = None
        try:
            response = self._session.get(header_url, timeout=20)
            if response.status_code == requests.codes.ok:
                headers = data_util.make_headers_from_string(response.text)
            else:
                headers = None
                self._logger.warning(
                    f'Error {response.status_code} when retrieving '
                    f'{header_url} headers.'
                )
        finally:
            if response is not None:
                response.close()
        self._logger.debug(f'End _retrieve_headers')
        return headers

    def _retrieve_json(self, source_name):
        # source name is a file id, because that's the only part that's
        # required to be unique for a retrieval from archive.gemini.edu
        #
        self._logger.debug(f'Begin _retrieve_json for {source_name}')
        gemini_url = f'{GEMINI_METADATA_URL}{source_name}'
        # Open the URL and fetch the JSON document for the observation
        response = None
        try:
            response = mc.query_endpoint_session(gemini_url, self._session)
            metadata = response.json()
        finally:
            if response is not None:
                response.close()
        if len(metadata) == 0:
            raise mc.CadcException(
                f'Could not find JSON record for {source_name} at '
                f'{gemini_url}.'
            )
        self._logger.debug(f'End _retrieve_json')
        return metadata

    def reset(self):
        super().reset()
        self._json_metadata = {}

    def set(self, storage_name):
        self.set_json_metadata(storage_name)
        self.set_headers(storage_name)

    def set_json_metadata(self, storage_name):
        """Retrieves Gemini JSON metadata to memory. Extracts a
        corresponding FileInfo record from that metadata."""
        for index, entry in enumerate(storage_name.destination_uris):
            if entry not in self._json_metadata.keys():
                self._logger.debug(f'Retrieve FileInfo for {entry}')
                temp = self._retrieve_json(storage_name.source_names[index])
                # the JSON that comes back is for multiple files, so find the
                # correct one
                for jj in temp:
                    # choose this key, and the comparison, because the lhs is
                    # a file id
                    if storage_name.source_names[index] in jj.get('filename'):
                        self._json_metadata[entry] = jj
                        f_name = self._json_metadata[entry].get('filename')
                        self._file_info[entry] = FileInfo(
                            id=entry,
                            size=self._json_metadata[entry].get('data_size'),
                            name=f_name,
                            md5sum=self._json_metadata[entry].get('file_md5'),
                            lastmod=self._json_metadata[entry].get('lastmod'),
                            file_type=data_util.get_file_type(f_name),
                            encoding=data_util.get_file_encoding(f_name),
                        )
                        break


class GeminiMetadataLookup:

    def __init__(self, metadata_reader):
        self._reader = metadata_reader

    def data_label(self, uri):
        temp = None
        if hasattr(self._reader, 'json_metadata'):
            temp = self._reader.json_metadata.get(uri).get('data_label')
        if temp is None:
            temp = self._reader.headers.get(uri)[0].get('DATALAB')
            if temp is None:
                temp = self._reader.headers.get(uri)[1].get('DATALAB')
        result = None
        if temp is not None:
            result = obs_file_relationship.repair_data_label(
                uri.split('/')[-1], temp
            )
        return result

    @property
    def reader(self):
        return self._reader

    @reader.setter
    def reader(self, value):
        self._reader = value
