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

import json
import logging
import os

from astropy.io.votable import parse_single_table
from astropy.table import Table
from bs4 import BeautifulSoup
from collections import OrderedDict
from datetime import datetime
from hashlib import md5

from caom2.diff import get_differences
from caom2pipe import manage_composable as mc

import gem2caom2.external_metadata as em
from gem2caom2 import composable


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
PLUGIN = os.path.join(os.path.dirname(THIS_DIR), 'main_app.py')
TEST_FILE = os.path.join(TEST_DATA_DIR, 'from_paul.txt')

FIRST_FILE_LIST = os.path.join(TEST_DATA_DIR, 'S20191010S0.html')
SECOND_FILE_LIST = os.path.join(TEST_DATA_DIR, 'S20191010S3.html')
STATE_FILE = '/usr/src/app/state.yml'
TOO_MANY_FILE_LIST = os.path.join(TEST_DATA_DIR, 'too_many_rows.html')

TEST_TODO_LIST = OrderedDict([(1572566500.951468, 'N20191101S0002.fits'),
                              (1572566505.241434, 'N20191101S0003.fits'),
                              (1572566511.547385, 'N20191101S0004.fits'),
                              (1572566515.857351, 'N20191101S0005.fits'),
                              (1572566666.679175, 'N20191101S0006.fits'),
                              (1572566670.978142, 'N20191101S0007.fits')])
TEST_BUILDER_OBS_ID = 'GS-2019B-Q-222-181-001'


# key - filename
# value - data label/observation id
TAP_QUERY_LOOKUP = {
    'S20161227S0051': 'GS-CAL20161227-5-001',
    'S20161227S0052': 'GS-CAL20161227-5-002',
    'S20161227S0053': 'GS-CAL20161227-5-003',
    'S20161227S0054': 'GS-CAL20161227-5-004',
    'S20161227S0055': 'GS-CAL20161227-5-005',
    'S20161227S0056': 'GS-CAL20161227-5-006',
    'S20161227S0057': 'GS-CAL20161227-5-007',
    '2004may20_0048': 'GS-CAL20040520-7-0048',
    '2004may20_0049': 'GS-CAL20040520-7-0049',
    '2004may20_0050': 'GS-CAL20040520-7-0050',
    '2004may20_0051': 'GS-CAL20040520-7-0051',
    '2004may20_0052': 'GS-CAL20040520-7-0052',
    '2004may20_0053': 'GS-CAL20040520-7-0053',
    '2004may20_0054': 'GS-CAL20040520-7-0054',
    '2004may20_0055': 'GS-CAL20040520-7-0055',
    '2004may20_0056': 'GS-CAL20040520-7-0056',
    '2004may20_0057': 'GS-CAL20040520-7-0057',
    '2004may20_0058': 'GS-CAL20040520-7-0058',
    '2004may20_0059': 'GS-CAL20040520-7-0059',
    '2004may20_0060': 'GS-CAL20040520-7-0060',
    '2004may20_0061': 'GS-CAL20040520-7-0061',
    '2004may20_0062': 'GS-CAL20040520-7-0062',
    'N20191101S0007': 'GN-2019B-ENG-1-160-008',
    'S20181219S0216': 'GS-CAL20181219-4-012',
    'S20181219S0217': 'GS-CAL20181219-4-013',
    'S20181219S0218': 'GS-CAL20181219-4-014',
    'S20181219S0219': 'GS-CAL20181219-4-015',
    'S20181219S0220': 'GS-CAL20181219-4-016',
    'S20181219S0221': 'GS-CAL20181219-4-017',
    'S20181219S0222': 'GS-CAL20181219-4-018',
    'S20181219S0223': 'GS-CAL20181219-4-019',
    'S20181219S0224': 'GS-CAL20181219-4-020',
    'S20181219S0225': 'GS-CAL20181219-4-021',
    'mrgN20060130S0149_add': 'GN-2006A-Q-90-1-001-MRG-ADD',
    'N20060130S0149': 'GN-2006A-Q-90-1-001',
    'N20060130S0150': 'GN-2006A-Q-90-1-002',
    'N20060130S0151': 'GN-2006A-Q-90-1-003',
    'N20200210S0077_bias': 'GN-CAL20200210-22-076-BIAS',
}


# structured by file id, observation id, filter_name (when looking up
# from SVOFPS for imaging files - x means there's no filter name), and
# instrument name
LOOKUP = {
    # bHROS
    'S20050825S0143': ['GS-2005B-SV-301-16-005', em.Inst.BHROS,
                       'GS-2005B-SV-301'],
    'S20051027S0089': ['GS-2005B-SV-302-20-001', em.Inst.BHROS,
                       'GS-2005B-SV-302'],
    'S20070130S0048': ['GS-2006B-Q-47-76-003', em.Inst.BHROS, 'GS-2006B-Q-47'],
    'S20070113S0060': ['GS-2006B-Q-7-32-008', em.Inst.BHROS, 'GS-2006B-Q-7'],
    # F2
    'S20150103S0144': ['GS-2014B-Q-17-53-030', em.Inst.F2, 'GS-2014B-Q-17'],
    'S20150103S0148': ['GS-2014B-Q-17-69-002', em.Inst.F2, 'GS-2014B-Q-17'],
    'S20150113S0186': ['GS-2014B-Q-60-11-012', em.Inst.F2, 'GS-2014B-Q-60'],
    'S20150101S0261': ['GS-2014B-Q-60-8-071', em.Inst.F2, 'GS-2014B-Q-60'],
    'S20150213S0110': ['GS-2015A-Q-60-126-001', em.Inst.F2, 'GS-2015A-Q-60'],
    'S20150622S0011': ['GS-2015A-Q-63-328-014', em.Inst.F2, 'GS-2015A-Q-63'],
    'S20150625S0077': ['GS-2015A-Q-63-365-016', em.Inst.F2, 'GS-2015A-Q-63'],
    'S20150625S0230': ['GS-2015A-Q-63-406-001', em.Inst.F2, 'GS-2015A-Q-63'],
    'S20160220S0125': ['GS-2016A-C-1-86-003', em.Inst.F2, 'GS-2016A-C-1'],
    'S20160616S0034': ['GS-2016A-FT-18-54-002', em.Inst.F2, 'GS-2016A-FT-18'],
    'S20171226S0358': ['GS-2017A-Q-28-280-012', em.Inst.F2, 'GS-2017A-Q-28'],
    'S20171123S0216': ['GS-2017B-Q-18-96-006', em.Inst.F2, 'GS-2017B-Q-18'],
    'S20171123S0166': ['GS-2017B-Q-45-156-079', em.Inst.F2, 'GS-2017B-Q-45'],
    'S20170221S0005': ['GS-CAL20170221-1-005', em.Inst.F2, 'GS-CAL20170221'],
    'S20181230S0026': ['GS-2018B-SV-301-144-020', em.Inst.F2,
                       'GS-2018B-SV-301'],
    'S20170905S0318': ['GS-2017A-Q-58-66-027', em.Inst.F2, 'GS-2017A-Q-58'],
    # Flamingos
    '02jul07.0186': ['GS-2002A-Q-13-2-0186', em.Inst.FLAMINGOS,
                     'GS-2002A-Q-13'],
    '02jun25.0071': ['GS-2002A-Q-7-1-0071', em.Inst.FLAMINGOS, 'GS-2002A-Q-7'],
    '02sep04.0230': ['GS-2002B-DD-2-6-0230', em.Inst.FLAMINGOS,
                     'GS-2002B-DD-2'],
    '02nov05.0072': ['GS-2002B-DS-1-31-0072', em.Inst.FLAMINGOS,
                     'GS-2002B-DS-1'],
    '02nov05.0011': ['GS-2002B-DS-1-41-0011', em.Inst.FLAMINGOS,
                     'GS-2002B-DS-1'],
    '02nov05.0016': ['GS-2002B-DS-1-41-0016', em.Inst.FLAMINGOS,
                     'GS-2002B-DS-1'],
    '02nov06.0031': ['GS-2002B-DS-1-46-0031', em.Inst.FLAMINGOS,
                     'GS-2002B-DS-1'],
    '02nov06.0042': ['GS-2002B-DS-1-46-0042', em.Inst.FLAMINGOS,
                     'GS-2002B-DS-1'],
    '02oct30.0205': ['GS-2002B-Q-5-20-0205', em.Inst.FLAMINGOS, 'GS-2002B-Q-5'],
    '02jun25.0115': ['GS-CAL20020625-10-0115', em.Inst.FLAMINGOS,
                     'GS-CAL20020625'],
    '02jun25.0121': ['GS-CAL20020625-11-0121', em.Inst.FLAMINGOS,
                     'GS-CAL20020625'],
    '02jul07.0034': ['GS-CAL20020707-8-0034', em.Inst.FLAMINGOS,
                     'GS-CAL20020707'],
    '02jul07.0040': ['GS-CAL20020707-9-0040', em.Inst.FLAMINGOS,
                     'GS-CAL20020707'],
    '02nov05.0001': ['GS-CAL20021105-1-0001', em.Inst.FLAMINGOS,
                     'GS-CAL20021105'],
    '02nov10.0008': ['GS-CAL20021110-2-0008', em.Inst.FLAMINGOS,
                     'GS-CAL20021110'],
    '02nov10.0174': ['GS-CAL20021110-6-0174', em.Inst.FLAMINGOS,
                     'GS-CAL20021110'],
    '02nov12.0455': ['GS-CAL20021112-26-0455', em.Inst.FLAMINGOS,
                     'GS-CAL20021112'],
    '02nov12.0024': ['GS-CAL20021112-3-0024', em.Inst.FLAMINGOS,
                     'GS-CAL20021112'],
    '02nov12.0043': ['GS-CAL20021112-4-0043', em.Inst.FLAMINGOS,
                     'GS-CAL20021112'],
    # GMOS
    'N20030107S0163': ['GN-2003A-Q-22-3-004', em.Inst.GMOS, 'GN-2003A-Q-22'],
    'N20071219S0193': ['GN-2007B-Q-112-14-018', em.Inst.GMOS, 'GN-2007B-Q-112'],
    'N20090313S0180': ['GN-2009A-Q-21-115-001', em.Inst.GMOS, 'GN-2009A-Q-21'],
    'N20100104S0208': ['GN-2009B-Q-121-15-001', em.Inst.GMOS, 'GN-2009B-Q-121'],
    'N20100104S0210': ['GN-2009B-Q-121-15-003', em.Inst.GMOS, 'GN-2009B-Q-121'],
    'N20100115S0346': ['GN-2010A-Q-35-10-002', em.Inst.GMOS, 'GN-2010A-Q-35'],
    'N20120105S0344': ['GN-2011A-Q-31-21-005', em.Inst.GMOS, 'GN-2011A-Q-31'],
    'N20131203S0006': ['GN-2013B-Q-28-150-002', em.Inst.GMOS, 'GN-2013B-Q-28'],
    'N20150217S0380': ['GN-2015A-C-2-96-002', em.Inst.GMOS, 'GN-2015A-C-2'],
    'N20150220S0320': ['GN-2015A-C-4-24-086', em.Inst.GMOS, 'GN-2015A-C-4'],
    'N20150216S0129': ['GN-2015A-Q-36-15-001', em.Inst.GMOS, 'GN-2015A-Q-36'],
    'N20150216S0142': ['GN-2015A-Q-91-5-002', em.Inst.GMOS, 'GN-2015A-Q-91'],
    'N20030104S0065': ['GN-CAL20030104-14-001', em.Inst.GMOS, 'GN-CAL20030104'],
    'N20030104S0161': ['GN-CAL20030104-18-003', em.Inst.GMOS, 'GN-CAL20030104'],
    'N20150217S0274': ['GN-CAL20150217-2-003', em.Inst.GMOS, 'GN-CAL20150217'],
    'N20150929S0013': ['GN-CAL20150925-2-007', em.Inst.GMOS, 'GN-CAL20150925'],
    'S20040124S0077': ['GS-2003B-Q-5-29-003', em.Inst.GMOS, 'GS-2003B-Q-5'],
    'S20060101S0075': ['GS-2005B-Q-22-17-006', em.Inst.GMOS, 'GS-2005B-Q-22'],
    'S20060103S0143': ['GS-2005B-Q-22-29-004', em.Inst.GMOS, 'GS-2005B-Q-22'],
    'S20060128S0316': ['GS-2005B-Q-27-33-001', em.Inst.GMOS, 'GS-2005B-Q-27'],
    'S20060122S0004': ['GS-2005B-Q-27-4-001', em.Inst.GMOS, 'GS-2005B-Q-27'],
    'S20060131S0110': ['GS-2005B-Q-54-2-004', em.Inst.GMOS, 'GS-2005B-Q-54'],
    'S20060620S0066': ['GS-2006A-Q-48-15-002', em.Inst.GMOS, 'GS-2006A-Q-48'],
    'S20080526S0024': ['GS-2008A-Q-41-26-013', em.Inst.GMOS, 'GS-2008A-Q-41'],
    'S20090620S0145': ['GS-2009A-Q-30-6-007', em.Inst.GMOS, 'GS-2009A-Q-30'],
    'S20060125S0027': ['GS-CAL20060125-1-002', em.Inst.GMOS, 'GS-CAL20060125'],
    'GN2001BQ013-04': ['GN2001BQ013-04', em.Inst.GMOS, 'GN-2001B-Q-13'],
    'N20110318S0581': ['GN-2011A-Q-87-69-001', em.Inst.GMOS, 'GN-2011A-Q-87'],
    # GNIRS
    'N20100915S0167': ['GN-2010B-Q-2-44-003', em.Inst.GNIRS, 'GN-2010B-Q-2'],
    'N20100722S0185': ['GN-2010B-SV-142-10-007', em.Inst.GNIRS,
                       'GN-2010B-SV-142'],
    'N20110323S0235': ['GN-2011A-Q-53-42-007', em.Inst.GNIRS, 'GN-2011A-Q-53'],
    'N20120104S0167': ['GN-2011B-Q-63-101-006', em.Inst.GNIRS, 'GN-2011B-Q-63'],
    'N20120102S0213': ['GN-2011B-Q-63-126-005', em.Inst.GNIRS, 'GN-2011B-Q-63'],
    'N20120101S0195': ['GN-2011B-Q-68-116-013', em.Inst.GNIRS, 'GN-2011B-Q-68'],
    'N20120103S0100': ['GN-2011B-Q-7-193-010', em.Inst.GNIRS, 'GN-2011B-Q-7'],
    'N20120103S0134': ['GN-2011B-Q-7-193-044', em.Inst.GNIRS, 'GN-2011B-Q-7'],
    'N20130419S0198': ['GN-2013A-Q-71-102-086', em.Inst.GNIRS, 'GN-2013A-Q-71'],
    'N20130408S0105': ['GN-2013A-Q-71-86-031', em.Inst.GNIRS, 'GN-2013A-Q-71'],
    'N20160123S0097': ['GN-2015B-SV-101-1061-005', em.Inst.GNIRS,
                       'GN-2015B-SV-101'],
    'N20151213S0022': ['GN-CAL20151213-6-002', em.Inst.GNIRS, 'GN-CAL20151213'],
    'N20160202S0098': ['GN-CAL20160202-3-039', em.Inst.GNIRS, 'GN-CAL20160202'],
    'S20041101S0215': ['GS-2004B-Q-19-20-023', em.Inst.GNIRS, 'GS-2004B-Q-19'],
    'N20170201S0246': ['GN-2017A-Q-44-25-031', em.Inst.GNIRS, 'GS-2017A-Q-44'],
    'N20180224S0063': ['GN-CAL20180224-3-001', em.Inst.GNIRS,
                       'GN-CAL20180224-3'],
    'N20171106S0187': ['GN-2017B-LP-16-470-002', em.Inst.GNIRS,
                       'GN-2017B-LP-16'],
    'N20170210S0013': ['GN-CAL20170209-5-003', em.Inst.GNIRS,
                       'GN-CAL20170209-5'],
    # GPI
    'S20140422S0167': ['GS-2014A-SV-408-6-003', em.Inst.GPI, 'GS-2014A-SV-408'],
    'S20180313S0108': ['GS-2018A-FT-101-5-043', em.Inst.GPI, 'GS-2018A-FT-101'],
    'S20140315S0348': ['GS-CAL20140315-4-004', em.Inst.GPI, 'GS-CAL20140315'],
    'S20140317S0028': ['GS-CAL20140316-5-012', em.Inst.GPI, 'GS-CAL20140316'],
    'S20140321S0011': ['GS-CAL20140320-7-011', em.Inst.GPI, 'GS-CAL20140320'],
    'S20140323S0255': ['GS-CAL20140323-6-014', em.Inst.GPI, 'GS-CAL20140323'],
    'S20150410S0541': ['GS-CAL20150410-1-007', em.Inst.GPI, 'GS-CAL20150410'],
    'S20160224S0323': ['GS-CAL20160224-12-017', em.Inst.GPI, 'GS-CAL20160224'],
    'S20150129S0028': ['GS-2014B-Q-501-10-001', em.Inst.GPI, 'GS-2014B-Q-501'],
    # GRACES
    'N20150604G0003': ['GN-CAL20150604-1000-1072', em.Inst.GRACES,
                       'GN-CAL20150604'],
    'N20150604G0014': ['GN-CAL20150604-1000-1081', em.Inst.GRACES,
                       'GN-CAL20150604'],
    'N20150807G0046': ['GN-CAL20150807-1000-1035', em.Inst.GRACES,
                       'GN-CAL20150807'],
    # GSAOI
    'S20130126S0134': ['GS-2012B-SV-499-21-002', em.Inst.GSAOI,
                       'GS-2012B-SV-499'],
    'S20140113S0002': ['GS-2013B-DD-1-13-002', em.Inst.GSAOI, 'GS-2013B-DD-1'],
    'S20140122S0227': ['GS-2013B-Q-26-19-004', em.Inst.GSAOI, 'GS-2013B-Q-26'],
    'S20140113S0167': ['GS-2013B-Q-61-8-008', em.Inst.GSAOI, 'GS-2013B-Q-61'],
    'S20130201S0246': ['GS-CAL20130201-3-017', em.Inst.GSAOI, 'GS-CAL20130201'],
    'S20140109S0210': ['GS-CAL20140109-3-009', em.Inst.GSAOI, 'GS-CAL20140109'],
    'S20181023S0087': ['GS-CAL20181023-5-001', em.Inst.GSAOI, 'GS-CAL20181023'],
    # HOKUPAA
    '2002APR23_591': ['GN-2002A-DD-1-319-591', em.Inst.HOKUPAA,
                      'GN-2002A-DD-1'],
    '2002APR24_007': ['GN-CAL20020424-1-007', em.Inst.HOKUPAA,
                      'GN-CAL20020424'],
    # hrwfs
    'S20030218S0027': ['GS-2003A-Q-6-1-002', em.Inst.HRWFS, 'GS-2003A-Q-6'],
    'S20030218S0042': ['GS-2003A-Q-6-1-002', em.Inst.HRWFS, 'GS-2003A-Q-6'],
    '2003jan05_0082': ['GS-CAL20030105-2-0072', em.Inst.HRWFS,
                       'GS-CAL20030105'],
    '2003mar03_0052': ['GS-CAL20030303-7-0006', em.Inst.HRWFS,
                       'GS-CAL20030303'],
    'S20030730S0036': ['GS-CAL20030730-10-006', em.Inst.HRWFS,
                       'GS-CAL20030730'],
    'S20031218S0049': ['GS-CAL20031218-1-034', em.Inst.HRWFS, 'GS-CAL20031218'],
    # Michelle
    'N20060705S0054': ['GN-2006A-C-14-49-002', em.Inst.MICHELLE,
                       'GN-2006A-C-14'],
    'N20060418S0123': ['GN-2006A-Q-58-9-007', em.Inst.MICHELLE,
                       'GN-2006A-Q-58'],
    'N20070310S0156': ['GN-2007A-C-11-234-001', em.Inst.MICHELLE,
                       'GN-2007A-C-11'],
    'N20080612S0038': ['GN-2008A-Q-43-6-002', em.Inst.MICHELLE,
                       'GN-2008A-Q-43'],
    'N20100131S0131': ['GN-2009B-C-1-62-001', em.Inst.MICHELLE, 'GN-2009B-C-1'],
    'N20110127S0219': ['GN-2010B-C-3-21-002', em.Inst.MICHELLE, 'GN-2010B-C-3'],
    'N20060413S0129': ['GN-CAL20060413-7-001', em.Inst.MICHELLE,
                       'GN-CAL20060413'],
    'N20080308S0086': ['GN-CAL20080308-8-016', em.Inst.MICHELLE,
                       'GN-CAL20080308'],
    # NICI
    'S20100102S0035': ['GS-2009B-Q-14-129-029', em.Inst.NICI, 'GS-2009B-Q-14'],
    'S20100218S0028': ['GS-2009B-Q-21-19-001', em.Inst.NICI, 'GS-2009B-Q-21'],
    'S20130216S0276': ['GS-2013A-Q-21-16-002', em.Inst.NICI, 'GS-2013A-Q-21'],
    'S20130528S0043': ['GS-2013A-Q-39-158-008', em.Inst.NICI, 'GS-2013A-Q-39'],
    'S20100110S0026': ['GS-CAL20100110-1-026', em.Inst.NICI, 'GS-CAL20100110'],
    'S20101227S0023': ['GS-2010B-Q-20-151-001', em.Inst.NICI, 'GS-2010B-Q-20'],
    # NIFS
    'N20120102S0663': ['GN-2011B-Q-59-70-016', em.Inst.NIFS, 'GN-2011B-Q-59'],
    'N20120411S0279': ['GN-2012A-Q-2-22-007', em.Inst.NIFS, 'GN-2012A-Q-2'],
    'N20120405S0322': ['GN-2012A-Q-57-151-005', em.Inst.NIFS, 'GN-2012A-Q-57'],
    'N20121229S0118': ['GN-2012B-Q-51-108-001', em.Inst.NIFS, 'GN-2012B-Q-51'],
    'N20120909S0132': ['GN-2012B-Q-73-172-002', em.Inst.NIFS, 'GN-2012B-Q-73'],
    'N20160426S0092': ['GN-2016A-Q-35-107-001', em.Inst.NIFS, 'GN-2016A-Q-35'],
    'N20161007S0382': ['GN-2016B-Q-27-31-007', em.Inst.NIFS, 'GN-2016B-Q-27'],
    'N20060723S0132': ['GN-2006A-C-11-670-006', em.Inst.NIFS, 'GN-2006A-C-11'],
    # NIRI
    'N20020405S0044': ['GN-2002A-C-4-2-001', em.Inst.NIRI, ''],
    'N20020620S0021': ['GN-2002A-C-5-1-001', em.Inst.NIRI, 'GN-2002A-C-5'],
    'N20020620S0035': ['GN-2002A-C-5-1-015', em.Inst.NIRI, 'GN-2002A-C-5'],
    'N20020620S0315': ['GN-2002A-C-5-21-002', em.Inst.NIRI, 'GN-2002A-C-5'],
    'N20150404S0726': ['GN-2015A-C-1-20-001', em.Inst.NIRI, 'GN-2015A-C-1'],
    'N20150404S0872': ['GN-2015A-C-1-27-001', em.Inst.NIRI, 'GN-2015A-C-1'],
    'N20150405S0028': ['GN-2015A-C-1-27-071', em.Inst.NIRI, 'GN-2015A-C-1'],
    'N20151129S0307': ['GN-2015B-Q-34-55-040', em.Inst.NIRI, 'GN-2015B-Q-34'],
    'N20090105S0057': ['GN-2008B-C-3-2-009', em.Inst.NIRI, 'GN-2008B-C-3'],
    'N20090105S0060': ['GN-2008B-C-3-2-012', em.Inst.NIRI, 'GN-2008B-C-3'],
    'N20090918S0386': ['GN-2009B-Q-52-332-011', em.Inst.NIRI, 'GN-20009B-Q-52'],
    'N20100620S0126': ['GN-2010A-Q-44-175-015', em.Inst.NIRI, 'GN-2009B-Q-44'],
    'N20120118S0441': ['GN-2011B-Q-28-10-002', em.Inst.NIRI, 'GN-2011B-Q-28'],
    'N20170302S0331': ['GN-2017A-SV-1-98-331', em.Inst.NIRI, 'GN-2017A-SV-1'],
    'N20170523S0331': ['GN-2017A-Q-70-214-001', em.Inst.NIRI, 'GN-2017A-Q-70'],
    'N20141203S0891': ['GN-2014B-C-1-157-100', em.Inst.NIRI, 'GN-2014B-C-1'],
    # OSCIR
    'r01dec05_007': ['GS-2001B-Q-31-9-007', em.Inst.OSCIR, 'GS-2001B-Q-31'],
    # Phoenix
    '2002jun10_0171': ['GS-2002A-DD-1-17-0171', em.Inst.PHOENIX,
                       'GS-2002A-DD-1'],
    '2003aug20_0073': ['GS-2003B-Q-51-27-0073', em.Inst.PHOENIX,
                       'GS-2003B-Q-51'],
    '2006apr07_0258': ['GS-2006A-C-10-1-0258', em.Inst.PHOENIX,
                       'GS-2006A-C-10'],
    '2006apr02_0073': ['GS-2006A-DD-1-1-0073', em.Inst.PHOENIX,
                       'GS-2006A-DD-1'],
    '2006dec10_0052': ['GS-2006B-C-8-2-0052', em.Inst.PHOENIX, 'GS-2006B-C-8'],
    '2002may12_0277': ['GS-CAL20020512-9-0277', em.Inst.PHOENIX,
                       'GS-CAL20020512'],
    # TReCS
    'S20050621S0037': ['GS-2005A-Q-15-1-001', em.Inst.TRECS, 'GS-2005A-Q-15'],
    'S20050918S0058': ['GS-2005B-Q-10-22-003', em.Inst.TRECS, 'GS-2005B-Q-10'],
    'S20080610S0045': ['GS-2008A-C-5-35-002', em.Inst.TRECS, 'GS-2008A-C-5'],
    'S20120922S0372': ['GS-2012A-Q-7-31-001', em.Inst.TRECS, 'GS-2012A-Q-7'],
    'S20050102S0024': ['GS-CAL20050102-1-001', em.Inst.TRECS, 'GS-CAL20050102'],
    'rS20050916S0159': ['GS-2003B-Q-23-17-001', em.Inst.TRECS,
                        'GS-2003B-Q-23'],
    # CIRPASS
    '2003mar09_1204': ['GS-2003A-Q-10-19-1204', em.Inst.CIRPASS,
                       'GS-2003A-Q-10'],
    '2003jun30_3507': ['GS-2003A-Q-14-2-3507', em.Inst.CIRPASS,
                       'GS-2003A-Q-14'],
    '2003mar09_1161': ['GS-2003A-Q-24-1-1161', em.Inst.CIRPASS,
                       'GS-2003A-Q-24'],
    '2003jul01_3598': ['GS-2003A-Q-3-22-3598', em.Inst.CIRPASS, 'GS-2003A-Q-3'],
    '2003mar08_1055': ['GS-CAL20030308-4-1055', em.Inst.CIRPASS,
                       'GS-CAL20030308'],
    '2003mar09_1247': ['GS-CAL20030309-9-1247', em.Inst.CIRPASS,
                       'GS-CAL20030309'],
    '2003mar13_1769': ['GS-CAL20030313-3-1769', em.Inst.CIRPASS,
                       'GS-CAL20030313'],
    '2003mar15_2027': ['GS-CAL20030315-18-2027', em.Inst.CIRPASS,
                       'GS-CAL20030315'],
    '2003mar20_2731': ['GS-CAL20030320-3-2731', em.Inst.CIRPASS,
                       'GS-CAL20030320'],
    '2003jun30_3384': ['GS-CAL20030630-1-3384', em.Inst.CIRPASS,
                       'GS-CAL20030630'],
    '2003jun30_3532': ['GS-CAL20030630-11-3532', em.Inst.CIRPASS,
                       'GS-CAL20030630'],
    '2003jun30_3463': ['GS-CAL20030630-4-3463', em.Inst.CIRPASS,
                       'GS-CAL20030630'],
    # TEXES
    'TX20071021_FLT.2037': ['GN-2007B-C-6-5-005-FLT', em.Inst.TEXES,
                            'GN-2007B-C-6'],
    'TX20131117_flt.3002': ['TX20131117_flt.3002', em.Inst.TEXES,
                            'GN-2013B-Q-38'],
    'TX20131117_raw.3002': ['TX20131117.3002', em.Inst.TEXES,
                            'GN-2013B-Q-38'],
    'TX20170321_flt.2505': ['TX20170321_flt.2505', em.Inst.TEXES,
                            'GN-2017A-Q-56'],
    'TX20170321_flt.2507': ['TX20170321_flt.2507', em.Inst.TEXES,
                            'GN-2017A-Q-56'],
    # processed
    'GS20141226S0203_BIAS': ['GS-CAL20141226-7-026-G-BIAS', em.Inst.GMOS,
                             'GS-CAL20141226'],
    'N20070819S0339_dark': ['GN-2007B-Q-107-150-004-DARK', em.Inst.GMOS,
                            'GN-2007B-Q-107'],
    'N20110927S0170_fringe': ['GN-CAL20110927-900-170-FRINGE', em.Inst.GMOS,
                              'GN-CAL20110927'],
    'N20120320S0328_stack_fringe': ['GN-CAL20120320-900-328-STACK-FRINGE',
                                    em.Inst.GMOS, 'GN-CAL20120320'],
    'N20130404S0512_flat': ['GN-2013A-Q-63-54-051-FLAT', em.Inst.NIRI,
                            'GN-2013A-Q-63'],
    'N20140313S0072_flat': ['GN-2013B-Q-75-163-011-FLAT', em.Inst.NIRI,
                            'GN-2013B-Q-75'],
    'N20141109S0266_bias': ['GN-CAL20141109-2-001-BIAS', em.Inst.GMOS,
                            'GN-CAL20141109'],
    'N20150804S0348_dark': ['GN-2015B-Q-53-138-061-DARK', em.Inst.GMOS,
                            'GN-2015B-Q-53'],
    'N20160403S0236_flat_pasted': ['GN-CAL20160404-7-017-FLAT-PASTED',
                                   em.Inst.GMOS, 'GN-CAL20160404'],
    'S20120922S0406': ['GS-2012B-Q-1-32-002', em.Inst.GMOS, 'GS-2012B-Q-1'],
    'S20131007S0067_fringe': ['GS-CAL20131007-900-067-FRINGE', em.Inst.GMOS,
                              'GS-CAL20131007'],
    'S20140124S0039_dark': ['GS-2013B-Q-16-277-019-DARK', em.Inst.F2,
                            'GS-2013B-Q-16'],
    'S20141129S0331_dark': ['GS-CAL20141129-1-001-DARK', em.Inst.F2,
                            'GS-CAL20141129'],
    'S20161227S0051': ['GS-CAL20161227-5-001', em.Inst.GMOS, 'GS-CAL20161227'],
    'fmrgN20020413S0120_add': ['GN-2002A-SV-78-7-003-FMRG-ADD', em.Inst.GMOS,
                               'GN-2002A-SV-78'],
    'gS20181219S0216_flat': ['GS-CAL20181219-4-021-G-FLAT', em.Inst.GMOS,
                             'GS-CAL20181219'],
    'gS20190301S0556_bias': ['GS-CAL20190301-4-046-G-BIAS', em.Inst.GMOS,
                             'GS-CAL20190301'],
    'mfrgS20041117S0073_add': ['GS-2004B-Q-42-1-001-MFRG-ADD', em.Inst.GMOS,
                               'GS-2004B-Q-42'],
    'mfrgS20160310S0154_add': ['GS-2016A-Q-7-175-001-MFRG-ADD', em.Inst.GMOS,
                               'GS-2016A-Q-7'],
    'mrgN20041016S0095': ['GN-2004B-Q-30-1-001', em.Inst.GMOS,
                          'GN-2004B-Q-30'],
    'mrgN20050831S0770_add': ['GN-2005B-Q-28-32-001-MRG-ADD', em.Inst.GMOS,
                              'GN-2005B-Q-28'],
    'mrgN20160311S0691_add': ['GN-2016A-Q-68-46-001-MRG-ADD', em.Inst.GMOS,
                              'GN-2016A-Q-68'],
    'mrgS20120922S0406': ['GS-2012B-Q-1-32-002', em.Inst.GMOS,
                          'GS-2012B-Q-1'],
    'mrgS20160901S0122_add': ['GS-2016B-Q-72-23-001-MRG-ADD', em.Inst.GMOS,
                              'GS-2016B-Q-72'],
    'mrgS20181016S0184_fringe': ['GS-CAL20181016-5-001-MRG-FRINGE',
                                 em.Inst.GMOS, 'GS-CAL20181016'],
    'rS20121030S0136': ['GS-2012B-Q-90-366-003', em.Inst.TRECS,
                        'GS-2012B-Q-90'],
    'rgS20100212S0301': ['GS-2010A-Q-36-5-246', em.Inst.GMOS,
                         'GS-2010A-Q-36'],
    'rgS20100316S0366': ['GS-2010A-Q-36-6-358', em.Inst.GMOS, 'GS-2010A-Q-36'],
    'rgS20130103S0098_FRINGE': ['GS-CAL20130103-3-001-RG-FRINGE', em.Inst.GMOS,
                                'GS-CAL20130103'],
    'rgS20131109S0166_FRINGE': ['GS-CAL20131109-17-001-RG-FRINGE',
                                em.Inst.GMOS, 'GS-CAL20131109'],
    'rgS20161227S0051_fringe': ['GS-CAL20161227-5-001-RG-FRINGE', em.Inst.GMOS,
                                'GS-CAL20161227'],
    'p2004may20_0048_FLAT': ['GS-CAL20040520-7-0048-P-FLAT', em.Inst.PHOENIX,
                             'GS-CAL20040520'],
    'p2004may19_0255_COMB': ['GS-2004A-Q-6-27-0255-P-COMB', em.Inst.PHOENIX,
                             'GS-2004A-Q-6'],
    'P2003JAN14_0148_DARK': ['GS-CAL20030114-7-0148', em.Inst.PHOENIX,
                             'GS-CAL2003011'],
    'P2002FEB03_0045_DARK10SEC': ['GS-CAL20020203-4-0045', em.Inst.PHOENIX,
                                  'GS-CAL20020203'],
    'P2002DEC02_0075_SUB.0001': ['GS-CAL20021202-3-0075', em.Inst.PHOENIX,
                                 'GS-CAL2002120'],
    '2004may19_0255': ['GS-2004A-Q-6-27-0255', em.Inst.PHOENIX, 'GS-2004A-Q-6'],
    'mrgN20060130S0149_add': ['GN-2006A-Q-90-1-001-MRG-ADD', em.Inst.GMOS,
                              'GS-2006A-Q-90'],
    'S20181016S0184': ['GS-CAL20181016-5-001', em.Inst.GMOS,
                       'GS-CAL20181016-5'],
    'N20200210S0077_bias': ['GN-CAL20200210-22-076-BIAS', em.Inst.GMOS,
                            'GN-CAL20200210'],
    'rgnN20140428S0171_flat': ['GN-2014A-Q-85-16-003-RGN-FLAT', em.Inst.NIFS,
                               'GN-2014A-Q-85']
}

call_count = 0


def mock_get_votable(url):
    try:
        x = url.split('/')
        filter_name = x[-1].replace('&VERB=0', '')
        votable = parse_single_table(
            f'{TEST_DATA_DIR}/votable/{filter_name}.xml')
        return votable, None
    except Exception as e:
        logging.error(f'get_vo_table failure for url {url}')
        logging.error(e)
        return None, None


def mock_get_pi_metadata(program_id):
    try:
        fname = f'{TEST_DATA_DIR}/programs/{program_id}.xml'
        with open(fname) as f:
            y = f.read()
            soup = BeautifulSoup(y, 'lxml')
            tds = soup.find_all('td')
            if len(tds) > 0:
                title = None
                if len(tds[1].contents) > 0:
                    title = tds[1].contents[0].replace('\n', ' ')
                pi_name = None
                if len(tds[3].contents) > 0:
                    pi_name = tds[3].contents[0]
                metadata = {'title': title,
                            'pi_name': pi_name}
                return metadata
        return None
    except Exception as e:
        logging.error(e)
        import traceback
        tb = traceback.format_exc()
        logging.error(tb)


def mock_get_file_info(archive, file_id):
    if '_prev' in file_id:
        return {'size': 10290,
                'md5sum': 'md5:{}'.format(
                    md5('-37'.encode()).hexdigest()),
                'type': 'image/jpeg',
                'name': file_id}
    else:
        return {'size': 665345,
                'md5sum': 'md5:a347f2754ff2fd4b6209e7566637efad',
                'type': 'application/fits',
                'name': file_id}


def mock_get_obs_metadata(file_id):
    try:
        fname = f'{TEST_DATA_DIR}/json/{file_id}.json'
        if os.path.exists(fname):
            with open(fname) as f:
                y = json.loads(f.read())
        else:
            # TODO
            y = [{'data_label': TAP_QUERY_LOOKUP.get(
                                    file_id, 'test_data_label'),
                  'filename': f'{file_id}.fits.bz2',
                  'lastmod': '2020-02-25T20:36:31.230'}]
        em.om.add(y, file_id)
    except Exception as e:
        logging.error(e)
        import traceback
        tb = traceback.format_exc()
        logging.error(tb)


class Object(object):
    pass

    def close(self):
        pass


def mock_query_endpoint(url, timeout=-1):
    # returns response.text
    result = Object()
    result.text = None
    global call_count

    if call_count == 0 and '20030106' not in url:
        with open(FIRST_FILE_LIST, 'r') as f:
            temp = f.read()
            now_dt = datetime.utcnow()
            now_date_str = datetime.strftime(now_dt, '%Y-%m-%d')
            now_time_str = datetime.strftime(now_dt, '%H:%M:%S')
            result.text = temp.replace(
                '2019-10-10', now_date_str).replace('05:09:24', now_time_str)
    elif call_count == 1 and '20030106' not in url:
        with open(SECOND_FILE_LIST, 'r') as f:
            result.text = f.read()
    elif '20030106' in url:
        with open(TOO_MANY_FILE_LIST, 'r') as f:
            result.text = f.read()
    else:
        if '20030107' in url or '20030105' in url:
            pass
        else:
            raise Exception('wut {} count {}'.format(url, call_count))
    call_count += 1
    return result


def mock_query_endpoint_3(url, timeout=-1):
    # returns json via response.text, depending on url
    result = Object()
    result.text = '[]'
    if 'jsonfilelist' in url:
        if 'filepre=S20200602' in url:
            with open(f'{TEST_DATA_DIR}/edu_query/S20200303_filepre.json',
                      'r') as f:
                result.text = f.read()
    return result


def mock_query_endpoint_2(url, timeout=-1):
    # returns response.json
    def x():
        if url.startswith('http://arcdev'):
            with open(f'{TEST_DATA_DIR}/jsonfilelist_composable_test.json',
                      'r') as f:
                temp = f.read()
        elif url == 'https://archive.gemini.edu/jsonsummary/canonical/' \
                    'notengineering/NotFail//filepre=N20191101S0001.fits':
            with open(f'{TEST_DATA_DIR}/json/N20191101S0001.json',
                      'r') as f:
                temp = f.read()
        else:
            # raise mc.CadcException(f'more wut? {url}')
            fid = url.split('filepre=')[1]
            # temp = [{"filename": "{fid}.bz2",
            #          "data_label": "GN-2019B-ENG-1-160-002",
            #          "lastmod": "2019-11-01 00:01:34.610517+00:00"}]
            temp = '[{"filename": "' + fid + '.bz2",' \
                   '"data_label": "GN-2019B-ENG-1-160-002",' \
                   '"lastmod": "2019-11-01 00:01:34.610517+00:00"}]'
            logging.error(f'wtf 2 {url}')
        return json.loads(temp)

    result = Object()
    result.json = x

    if url.startswith('http://arcdev'):
        pass
    elif url == 'https://archive.gemini.edu/jsonsummary/canonical/' \
                'notengineering/NotFail//filepre=N20191101S0001.fits':
        pass
    else:
        # raise mc.CadcException(f'wut? {url}')
        pass
    return result


def mock_write_state(start_time):
    test_bookmark = {'bookmarks': {composable.GEM_BOOKMARK:
                                   {'last_record': start_time}}}
    mc.write_as_yaml(test_bookmark, STATE_FILE)


def mock_write_state2(prior_timestamp=None):
    # to ensure at least one spin through the execution loop, test case
    # must have a starting time greater than one config.interval prior
    # to 'now', default interval is 10 minutes
    if prior_timestamp is None:
        prior_s = datetime.utcnow().timestamp() - 15 * 60
    else:
        prior_s = mc.make_seconds(prior_timestamp)
    test_start_time = datetime.fromtimestamp(prior_s)
    test_bookmark = {'bookmarks': {composable.GEM_BOOKMARK:
                                   {'last_record': test_start_time}}}
    mc.write_as_yaml(test_bookmark, STATE_FILE)


def mock_repo_create(arg1):
    # arg1 is an Observation instance
    act_fqn = f'{TEST_DATA_DIR}/{arg1.observation_id}.actual.xml'
    ex_fqn = f'{TEST_DATA_DIR}/{arg1.observation_id}.expected.xml'
    mc.write_obs_to_file(arg1, act_fqn)
    result = compare(ex_fqn, act_fqn, arg1.observation_id)
    if result is not None:
        assert False, result


read_call_count = 0


def mock_repo_read(arg1, arg2):
    logging.error(f'arg1 {arg1} arg2 {arg2}')
    # arg1 GEMINI arg2 GS-CAL20191010-3-034
    global read_call_count
    if arg1 == 'GEMINI' and arg2 == 'GS-2004A-Q-6-27-0255':
        return ''
    elif read_call_count == 0:
        read_call_count = 1
        return None
    else:
        return mc.read_obs_from_file(
            f'{TEST_DATA_DIR}/GS-2019B-Q-222-181-001.expected.xml')


def mock_repo_update(ignore1):
    return None


def compare(ex_fqn, act_fqn, entry):
    ex = mc.read_obs_from_file(ex_fqn)
    act = mc.read_obs_from_file(act_fqn)
    result = get_differences(ex, act, 'Observation')
    if result:
        result_str = '\n'.join([r for r in result])
        msg = f'Differences found obs id {ex.observation_id} ' \
              f'file id {entry} ' \
              f'instr {ex.instrument.name}\n{result_str}'
        return msg
    return None


def _query_mock_none(ignore1, ignore2):
    return Table.read('observationID,lastModified\n'.split('\n'), format='csv')


def _query_mock_one(ignore1, ignore2):
    return Table.read('observationID,lastModified\n'
                      'test_data_label,2020-02-25T20:36:31.230\n'.split('\n'),
                      format='csv')


def mock_query_tap(query_string, mock_tap_client):
    if query_string.startswith('SELECT A.uri'):
        return Table.read(
            f'uri,lastModified\n'
            f'gemini:GEMINI/N20191101S0007.fits,'
            f'2020-02-25T20:36:31.230\n'.split('\n'), format='csv')
    else:
        file_id = query_string.split(
            'GEMINI/')[1].replace('\'', '').replace('.fits', '').strip()
        result = TAP_QUERY_LOOKUP.get(file_id, 'test_data_label')
        return Table.read(f'observationID,lastModified\n'
                          f'{result},2020-02-25T20:36:31.230\n'.split('\n'),
                          format='csv')
