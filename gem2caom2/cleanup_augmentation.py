# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
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

import logging

from caom2 import Observation
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from gem2caom2 import obs_file_relationship


def visit(observation, **kwargs):
    """
    If there are artifacts with the same name, but different case, prefer
    the lower case artifact, and remove the upper-case one.

    :param observation: Observation instance - check all it's artifacts
    :param kwargs:
    """
    mc.check_param(observation, Observation)
    artifact_count = 0
    plane_count = 0

    if len(observation.planes.values()) > 1:
        all_artifact_keys = cc.get_all_artifact_keys(observation)
        all_artifact_keys_lower = [ii.lower() for ii in all_artifact_keys]
        set_artifact_keys_lower = set(all_artifact_keys_lower)
        delete_these_artifacts = []
        if len(all_artifact_keys) != len(set_artifact_keys_lower):
            for entry in set_artifact_keys_lower:
                ignore_scheme, ignore_path, file_name = mc.decompose_uri(entry)
                file_id = obs_file_relationship.remove_extensions(file_name)
                # it's the suffix that has the different case, so use it
                # to figure out which artifacts shouldn't exist
                suffixes = obs_file_relationship.get_suffix(
                    file_id, observation.observation_id
                )
                for key in all_artifact_keys:
                    for suffix in suffixes:
                        if suffix.upper() in key:
                            # get the fits, previews, thumbnails as well
                            delete_these_artifacts.append(key)

        delete_these_planes = []
        for entry in delete_these_artifacts:
            for plane in observation.planes.values():
                if entry in plane.artifacts.keys():
                    plane.artifacts.pop(entry)
                    logging.info(f'Removing {entry} from {plane.product_id}.')
                    artifact_count += 1
                if len(plane.artifacts.keys()) == 0:
                    delete_these_planes.append(plane.product_id)

        for entry in set(delete_these_planes):
            observation.planes.pop(entry)
            logging.info(
                f'Removing {entry} from {observation.observation_id}.'
            )
            plane_count += 1

    logging.debug(
        f'Completed cleanup for {observation.observation_id}. Removed '
        f'{artifact_count} artifacts and {plane_count} planes.'
    )
    result = {
        'artifacts': artifact_count,
        'planes': plane_count,
    }
    return observation
