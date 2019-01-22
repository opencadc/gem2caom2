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

import io
import logging
import re

import requests
from astropy.io.votable import parse_single_table
from caom2pipe import manage_composable as mc
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

SVO_URL = 'http://svo2.cab.inta-csic.es/svo/theory/fps/fps.php?ID=Gemini/'


def get_votable(url):
    """
    Download the VOTable XML for the given url and return a astropy.io.votable
    object.

    :param url: query url for the SVO service
    :return: astropy.io.votable of the first table and
             an error_message if there was an error downloading the data
    """
    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    votable = None
    response = None
    error_message = None
    try:
        response = session.get(url)
        fh = io.BytesIO(bytes(response.text, 'utf-8'))
        response.close()
        votable = parse_single_table(fh)
    except Exception as e:
        error_message = str(e)
    if response:
        response.close()
    return votable, error_message


def filter_metadata(instrument, filters):
    """
    For the given instrument and filters, go to the SVO Filter Profile Service

    http://svo2.cab.inta-csic.es/svo/theory/fps/

    and return energy metadata.

    :param instrument: The instrument name.
    :param filters: The filter name.
    :return: Energy metadata dictionary.
    """

    try:
        filter_md = {}
        filter_names = filters.split('&')
        wl_min = 0.0
        wl_max = 100000.0
        width_min = 100000.0

        for filter_name in filter_names:
            if 'Hartmann' in filter_name:
                continue
            if filter_name == 'open':
                if 'GMOS' in instrument:
                    wmin = 3500.0
                    wmax = 11000.0
                    wl_width = wmax - wmin
                    wl_eff = (wmax + wmin)/2.0
            elif filter_name == 'GG455':
                wmin = 4600.0
                wmax = 11000.0
                wl_width = wmax - wmin
                wl_eff = (wmax + wmin)/2.0
            elif filter_name == 'OG515':
                wmin = 5200.0
                wmax = 11000.0
                wl_width = wmax - wmin
                wl_eff = (wmax + wmin)/2.0
            elif filter_name == 'RG610':
                wmin = 6150.0
                wmax = 11000.0
                wl_width = wmax - wmin
                wl_eff = (wmax + wmin)/2.0
            elif filter_name == 'RG780':
                wmin = 780.0
                wmax = 11000.0
                wl_width = wmax - wmin
                wl_eff = (wmax + wmin)/2.0
            else:
                if instrument == 'F2':
                    instrument = 'Flamingos2'
                if instrument == 'NIRI':
                    filter_name = re.sub(r'con', 'cont', filter_name)
                    # SVO filter service has renamed some Gemini NIRI filters...
                    if filter_name == 'H2v=2-1s1-G0220':
                        filter_name = 'H2S1v2-1-G0220'
                    if filter_name == 'H2v=1-0s1-G0216':
                        filter_name = 'H2S1v1-0-G0216'

                filter_id = "{}.{}".format(instrument, filter_name)
                url = "{}{}".format(SVO_URL, filter_id)

                # Open the URL and fetch the VOTable document.
                # Some Gemini filters in SVO filter database have bandpass info
                # only for 'w'arm filters.  First check for filter without 'w'
                # appended to the ID (which I assume means bandpass is for cold
                # filter), then search for 'w' if nothing is found...
                logging.error(filter_id)
                votable, error_message = get_votable(url)
                if not votable:
                    url += 'w'
                    votable, error_message = get_votable(url)
                if not votable:
                    raise mc.CadcException(
                        'Unable to download SVO filter data from {} because {}'
                            .format(url, error_message))

                wl_width = votable.get_field_by_id('WidthEff').value
                wl_eff = votable.get_field_by_id('WavelengthEff').value
                wmin = wl_eff - wl_width/2.0
                wmax = wl_eff + wl_width/2.0

            if wmin > wl_min:
                wl_min = wmin
            if wmax < wl_max:
                wl_max = wmax
            if wl_width < width_min:
                width_min = wl_width

        filter_md['wl_min'] = wmin
        filter_md['wl_max'] = wmax
        filter_md['wl_eff_width'] = wl_width
        filter_md['wl_eff'] = wl_eff
        logging.debug('\tFilter(s): {}  MD: {}'.format(filter_names, filter_md))
        return filter_md
    except Exception as e:
        logging.error(e)
