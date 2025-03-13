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
#  : 4 $
#
# ***********************************************************************
#

import logging

from bs4 import BeautifulSoup

from caom2pipe import manage_composable as mc


__all__ = ['MDContext', 'PIMetadata']


class PIMetadata:
    """Store the PI Metadata query results for the length of the application run."""

    def __init__(self, gemini_session):
        self._pm = {}
        self._gemini_session = gemini_session
        self._logger = logging.getLogger(self.__class__.__name__)

    def get_pi_metadata(self, program_id):
        self._logger.debug(f'Begin get_pi_metadata for {program_id}')
        metadata = self._pm.get(program_id)
        if not metadata:
            # for TaskType.SCRAPE
            if self._gemini_session is None:
                metadata = None
                self._logger.warning(f'No external access. No PI metadata.')
            else:
                program_url = f'https://archive.gemini.edu/programinfo/{program_id}'
                # Open the URL and fetch the JSON document for the observation
                response = None
                try:
                    self._logger.debug(f'Querying {program_url} for program metadata.')
                    response = mc.query_endpoint_session(program_url, self._gemini_session)
                    xml_metadata = response.text
                finally:
                    if response:
                        response.close()
                metadata = {}
                soup = BeautifulSoup(xml_metadata, 'lxml')
                tds = soup.find_all('td')
                if len(tds) > 0:
                    # sometimes the program id points to an html page with an
                    # empty table, see e.g. N20200210S0077_bias
                    title = None
                    if len(tds[1].contents) > 0:
                        title = tds[1].contents[0].replace('\n', ' ')
                    pi_name = None
                    if len(tds[3].contents) > 0:
                        pi_name = tds[3].contents[0]
                    metadata = {
                        'title': title,
                        'pi_name': pi_name,
                    }
                # keep empty entries to minimize archive.gemini.edu queries each run
                self._pm[program_id] = metadata
        self._logger.debug('End get_pi_metadata')
        return metadata


class MDContext:
    """Composition class to hold the filter cache and pi metadata for easier attachment to GemName instances."""

    def __init__(self, filter_cache, pi_metadata):
        self.filter_cache = filter_cache
        self.pi_metadata = pi_metadata
