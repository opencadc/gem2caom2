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

from caom2pipe import manage_composable as mc
from caom2pipe import execute_composable as ec
from gem2caom2.gem_name import GemName


class GeminiObsMetadata(object):
    """A place to hold access to output from multiple jsonsummary
    queries.

    Hold the query results for all files associated with an observation.
    Use the 'add' method to add a single jsonsummary query result.
    Use the 'reset_index' method to have the 'get' method look up the
    results associated with a particular file_id.
    """

    def __init__(self):
        # a dictionary of all the jsonsummary results
        self.lookup = {}
        # which dictionary entry is of current lookup interest
        self.current = None
        # the json summary results are a list, track which entry in the
        # list has the information for a particular file_id
        self.index = -1

    def add(self, metadata, file_id):
        self.lookup[file_id] = metadata
        self._reset_index(file_id)

    def get(self, lookup):
        return self.current[self.index].get(lookup)

    def reset_index(self, uri):
        file_id = GemName.remove_extensions(ec.CaomName(uri).file_name)
        self._reset_index(file_id)

    def _reset_index(self, file_id):
        if file_id not in self.lookup:
            raise mc.CadcException(
                'ObsMetadata: Mystery file id {}'.format(file_id))
        self.current = self.lookup[file_id]
        self.index = self._get_index(file_id)

    def _get_index(self, file_id):
        result = -1
        for index, value in enumerate(self.current):
            indexed_f_name = value.get('filename')
            if indexed_f_name is not None:
                temp = GemName.remove_extensions(indexed_f_name)
                if temp == file_id:
                    result = index
                    break
        if result == -1:
            # TODO - set obs id?
            raise mc.CadcException(
                'JSON Summary: unrecognized file_id {} in obs_id {}'.format(
                    file_id, ''))
        return result
