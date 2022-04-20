# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
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

from datetime import datetime

from caom2 import Observation
from caom2pipe import manage_composable as mc

FILE_URL = 'https://archive.gemini.edu/file'


def visit(observation, **kwargs):
    """
    If the observation says the data release date is past, attempt to
    retrieve the fits file if it is not already at CADC.
    """
    mc.check_param(observation, Observation)
    working_dir = kwargs.get('working_directory', './')
    clients = kwargs.get('clients')
    if clients is None:
        logging.warning('Need clients to update. Stopping pull visitor.')
        return
    observable = kwargs.get('observable')
    if observable is None:
        raise mc.CadcException('Visitor needs a observable parameter.')
    metadata_reader = kwargs.get('metadata_reader')
    if metadata_reader is None:
        raise mc.CadcException('Visitor needs a metadata_reader parameter.')
    storage_name = kwargs.get('storage_name')
    if storage_name is None:
        raise mc.CadcException('Visitor needs a storage_name parameter.')

    count = 0
    if observable.rejected.is_bad_metadata(observation.observation_id):
        logging.info(
            f'Stopping visit for {observation.observation_id} '
            f'because of bad metadata.'
        )
    else:
        for plane in observation.planes.values():
            if (
                plane.data_release is None
                or plane.data_release > datetime.utcnow()
            ):
                logging.info(
                    f'Plane {plane.product_id} is proprietary. No file '
                    f'access.'
                )
                continue

            for artifact in plane.artifacts.values():
                # compare file names, because part of this visitor is to
                # change the URIs
                artifact_f_name = artifact.uri.split('/')[-1]
                if artifact_f_name != storage_name.file_name:
                    logging.debug(
                        f'Leave {artifact.uri}, want {storage_name.file_uri}'
                    )
                    continue
                try:
                    f_name = mc.CaomName(artifact.uri).file_name
                    if '.jpg' not in f_name:
                        logging.debug(f'Checking for {f_name}')
                        file_url = f'{FILE_URL}/{f_name}'
                        fqn = os.path.join(working_dir, f_name)

                        # want to compare the checksum from the JSON, and the
                        # checksum at CADC storage - if they are not the same,
                        # retrieve the file from archive.gemini.edu again
                        json_md5sum = metadata_reader.file_info.get(
                            artifact.uri
                        ).md5sum
                        look_pull_and_put(
                            artifact.uri, fqn, file_url, clients, json_md5sum
                        )
                        if os.path.exists(fqn):
                            logging.info(
                                f'Removing local copy of {f_name} after '
                                f'successful storage call.'
                            )
                            os.unlink(fqn)
                except Exception as e:
                    if not (
                        observable.rejected.check_and_record(
                            str(e), observation.observation_id
                        )
                    ):
                        raise e
    logging.info(f'Completed pull visitor for {observation.observation_id}.')
    result = {'observation': count}
    return observation


def look_pull_and_put(storage_name, fqn, url, clients, checksum):
    """Checks to see if a file exists at CADC. If yes, stop. If no,
    pull via https to local storage, then put to CADC storage.

    :param storage_name Artifact URI as the file will appear at CADC
    :param fqn name on disk for caching between the
        pull and the put
    :param url for retrieving the file externally, if it does not exist
    :param clients GemClientCollection instance
    :param checksum what the CAOM observation says the checksum should be -
        just the checksum part of ChecksumURI please, or the comparison will
        always fail.
    """
    cadc_meta = clients.data_client.info(storage_name)
    if (
        checksum is not None
        and cadc_meta is not None
        and cadc_meta.md5sum.replace('md5:', '') != checksum
    ) or cadc_meta is None:
        logging.debug(
            f'Different checksums: Source {checksum}, CADC {cadc_meta}'
        )
        mc.http_get(url, fqn)
        clients.data_client.put(os.path.dirname(fqn), storage_name)
        logging.info(
            f'Retrieved {os.path.basename(fqn)} for storage as '
            f'{storage_name}'
        )
    else:
        logging.info(f'{os.path.basename(fqn)} already exists at CADC.')
