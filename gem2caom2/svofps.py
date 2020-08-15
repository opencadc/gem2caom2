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

from caom2pipe import astro_composable as ac


def filter_metadata(instrument, filters):
    """
    For the given instrument and filters, go to the SVO Filter Profile Service

    http://svo2.cab.inta-csic.es/svo/theory/fps/

    and return energy metadata.

    :param instrument: The instrument name.
    :param filters: The filter name.
    :return: FilterMetadata instance, or None, if there's no SVO information
        for the filter.
    """

    try:
        filter_names = filters.split('+')
        # use detector maximums as defaults
        wl_min = 0.0
        wl_max = 100000.0
        width_min = 100000.0
        wl_width = wl_max - wl_min
        wl_eff = (wl_max + wl_min)/2.0

        # does the filter exist at SVO?
        filter_name_found = False

        for index in filter_names:
            filter_name = index.strip()
            filter_id = "{}.{}".format(instrument, filter_name)
            # VERB=0 parameter means the smallest amount returned
            if instrument == 'Flamingos':
                url = "{}KPNO/{}&VERB=0".format(ac.SVO_URL, filter_id)
            else:
                url = "{}Gemini/{}&VERB=0".format(ac.SVO_URL, filter_id)

            # Open the URL and fetch the VOTable document.
            # Some Gemini filters in SVO filter database have bandpass info
            # only for 'w'arm filters.  First check for filter without 'w'
            # appended to the ID (which I assume means bandpass is for cold
            # filter), then search for 'w' if nothing is found...
            votable, error_message = ac.get_vo_table(url)
            if not votable:
                if instrument == 'Flamingos':
                    url = "{}KPNO/{}w&VERB=0".format(ac.SVO_URL, filter_id)
                else:
                    url = "{}Gemini/{}w&VERB=0".format(ac.SVO_URL, filter_id)
                votable, error_message = ac.get_vo_table(url)
            if not votable:
                logging.error(
                    'Unable to download SVO filter data from {} because {}'
                    .format(url, error_message))
                continue

            # DB - 14-04-19 After discussion with a few others use the
            # wavelength lookup values “WavelengthCen” and “FWHM” returned
            # from the SVO. Looking at some of the IR filters
            # use the more common “WavelengthCen” and “FWHM” values that
            # the service offers.

            filter_name_found = True
            wl_width = votable.get_field_by_id('FWHM').value
            wl_eff = votable.get_field_by_id('WavelengthCen').value
            w_min = wl_eff - wl_width/2.0
            w_max = wl_eff + wl_width/2.0

            if w_min > wl_min:
                wl_min = w_min
            if w_max < wl_max:
                wl_max = w_max
            if wl_width < width_min:
                width_min = wl_width

        if filter_name_found:
            fm = FilterMetadata(instrument)
            # SVO filter units are angstroms, Gemini CAOM2 spectral wcs is
            # microns
            fm.central_wl = wl_eff / 1.0e4
            fm.bandpass = wl_width / 1.0e4
            logging.info(
                'Filter(s): {}  MD: {}, {}'.format(filter_names, fm.central_wl,
                                                   fm.bandpass))
            return fm
        else:
            return None
    except Exception as e:
        logging.error(e)
        import traceback
        tb = traceback.format_exc()
        logging.error(tb)


class FilterMetadata(object):
    # DB - 27-02-19
    # I’ve approached Gemini spectra like I have for DAO but the header
    # content for DAO spectra permits an accurate estimate of the upper/lower
    # wavelengths (a few % accuracy, at least for data from the last decade
    # or so).  It’s more a ‘ballpark’ estimate for Gemini even for simple
    # long-slit instruments like PHOENIX and NIRI.  And this approach isn’t
    # very consistent for spectrographs like GRACES or bHROS or GPI.  e.g.
    # for GRACES the WCS I’ve built  suggests that the full spectral range
    # is covered from one end of the detector to the other but in reality
    # there are perhaps 40 individual spectral orders stacked one above
    # the other each covering a few % of the full wavelength range with a
    # bit of overlap between successive orders.  I’ve never actually tried
    # to do a spectral cutout but the returned cutouts for unprocessed
    # GRACES data wouldn’t be valid.
    #
    # Just talked with Chris briefly as well.  For consistency and to avoid
    # misleading users it’s likely best to use the ‘range’ approach.
    #
    # microns for the units, not nm.  And second line would be something like:
    #
    # axis.range = CoordRange1D(RefCoord(0.5,crval1),
    #                           RefCoord(1.5,(crval1+bandpass)))
    #
    # and crval1 is the lower wavelength.

    def __init__(self, instrument=None):
        self.central_wl = None
        self.bandpass = None
        self.resolving_power = None
        self.instrument = instrument

    @property
    def central_wl(self):
        """Central wavelength for a filter."""
        return self._central_wl

    @central_wl.setter
    def central_wl(self, value):
        self._central_wl = value

    @property
    def bandpass(self):
        """Width of a filter."""
        return self._bandpass

    @bandpass.setter
    def bandpass(self, value):
        self._bandpass = value

    @property
    def resolving_power(self):
        if self.instrument in ['Phoenix', 'TEXES']:
            return None
        elif self._resolving_power is None:
            if self.instrument == 'NIRI':
                return None
            else:
                self.adjust_resolving_power()
                return self._resolving_power
        else:
            return self._resolving_power

    @resolving_power.setter
    def resolving_power(self, value):
        self._resolving_power = value

    def adjust_bandpass(self, variance):
        self.bandpass = ((self.central_wl + variance) -
                         (self.central_wl - variance))

    def set_bandpass(self, w_max, w_min):
        self.bandpass = (w_max - w_min)

    def set_central_wl(self, w_max, w_min):
        self.central_wl = (w_max + w_min) / 2.0

    def set_resolving_power(self, w_max, w_min):
        self.resolving_power = (w_max + w_min) / (2 * self.bandpass)

    def adjust_resolving_power(self):
        # the formula for direct imaging data
        self.resolving_power = self.central_wl / self.bandpass
