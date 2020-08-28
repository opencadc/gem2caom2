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

from datetime import datetime

import matplotlib.image as image

from caom2 import Observation, ProductType, ReleaseType
from caom2pipe import manage_composable as mc
from gem2caom2.gem_name import GemName, ARCHIVE

__all__ = ['visit']


PREVIEW_URL = 'https://archive.gemini.edu/preview/'
MIME_TYPE = 'image/jpeg'


def visit(observation, **kwargs):
    mc.check_param(observation, Observation)

    working_dir = kwargs.get('working_directory', './')
    cadc_client = kwargs.get('cadc_client')
    if cadc_client is None:
        logging.warning('Need a cadc_client to update preview records.')
    stream = kwargs.get('stream')
    if stream is None:
        raise mc.CadcException('Visitor needs a stream parameter.')
    observable = kwargs.get('observable')
    if observable is None:
        raise mc.CadcException('Visitor needs a observable parameter.')

    count = 0
    for plane in observation.planes.values():
        if (plane.data_release is None or
                plane.data_release > datetime.utcnow()):
            logging.info(f'Plane {plane.product_id} is proprietary. No '
                         f'preview access or thumbnail creation.')
            continue
        count += _do_prev(observation.observation_id, working_dir,
                          plane, cadc_client, stream, observable)
    logging.info('Completed preview augmentation for {}.'.format(
        observation.observation_id))
    return {'artifacts': count}


def _check_for_delete(file_name, uri, observable, plane):
    """If the preview file doesn't exist, but the artifact that represents it
    does, remove that artifact from the Observation instance."""
    result = 0
    if (observable.rejected.is_no_preview(
            file_name) and uri in plane.artifacts.keys()):
        logging.warning(f'Removing artifact for non-existent preview {uri}')
        plane.artifacts.pop(uri)
        result = 1
    return result


def _do_prev(obs_id, working_dir, plane, cadc_client, stream, observable):
    """Retrieve the preview file, so that a thumbnail can be made,
    store the preview if necessary, and the thumbnail, to ad.
    Then augment the CAOM observation with the two additional artifacts.
    """
    count = 0
    gem_name = GemName(obs_id=obs_id, file_id=plane.product_id)
    preview = gem_name.prev
    if observable.rejected.is_no_preview(preview):
        logging.info(f'Stopping visit because no preview exists for {preview} '
                     f'in observation {obs_id}.')
        observable.rejected.record(mc.Rejected.NO_PREVIEW, preview)
        count += _check_for_delete(
            preview, gem_name.prev_uri, observable, plane)
    else:
        preview_fqn = os.path.join(working_dir, preview)
        thumb = gem_name.thumb
        thumb_fqn = os.path.join(working_dir, thumb)

        # get the file - try disk first, then CADC, then Gemini
        if not os.access(preview_fqn, 0) and cadc_client is not None:
            try:
                mc.data_get(cadc_client, working_dir, preview, ARCHIVE,
                            observable.metrics)
            except mc.CadcException:
                preview_url = '{}{}.fits'.format(PREVIEW_URL, plane.product_id)
                try:
                    mc.http_get(preview_url, preview_fqn)
                    mc.data_put(cadc_client, working_dir, preview, ARCHIVE,
                                stream, MIME_TYPE, mime_encoding=None,
                                metrics=observable.metrics)
                except Exception as e:
                    if observable.rejected.check_and_record(str(e), preview):
                        count += _check_for_delete(
                            preview, gem_name.prev_uri, observable, plane)
                        return count
                    else:
                        raise e

        try:
            fp = open(preview_fqn, 'r')
        except PermissionError as e:
            raise mc.CadcException(
                f'Should not have reached this point in thumbnail generation '
                f'for {plane.product_id}')

        _augment(plane, gem_name.prev_uri, preview_fqn, ProductType.PREVIEW)
        count = 1

        logging.debug(f'Generate thumbnail for file id {plane.product_id}')
        if os.access(thumb_fqn, 0):
            os.remove(thumb_fqn)
        thumb_fig = image.thumbnail(preview_fqn, thumb_fqn, scale=0.25)
        if thumb_fig is not None:
            count = 1

        thumb_uri = gem_name.thumb_uri
        _augment(plane, thumb_uri, thumb_fqn, ProductType.THUMBNAIL)
        if cadc_client is not None:
            mc.data_put(cadc_client, working_dir, thumb, ARCHIVE, stream,
                        MIME_TYPE, mime_encoding=None,
                        metrics=observable.metrics)
        count += 1
    return count


def _augment(plane, uri, fqn, product_type):
    temp = None
    if uri in plane.artifacts:
        temp = plane.artifacts[uri]
    plane.artifacts[uri] = mc.get_artifact_metadata(
        fqn, product_type, ReleaseType.DATA, uri, temp)
