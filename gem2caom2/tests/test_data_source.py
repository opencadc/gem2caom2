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

from datetime import datetime
from mock import Mock, patch
from gem2caom2 import data_source
import gem_mocks


@patch('caom2pipe.manage_composable.query_endpoint_session')
def test_incremental_source(query_mock):
    # clean up the singleton cache content as initial conditions
    json_lookup.flush()

    # https://archive.gemini.edu/jsonsummary/canonical/entrytimedaterange=
    # 2021-01-01T20:03:00.000000%202021-01-01T22:13:00.000000/
    # ?orderby=entrytime
    # get results
    query_mock.side_effect = gem_mocks.mock_query_endpoint_2

    test_subject = data_source.IncrementalSource(Mock())
    assert len(test_subject._json_cache._lookup) == 0, 'initial conditions'
    assert test_subject is not None, 'expect construction success'
    prev_exec_time = datetime(
        year=2021, month=1, day=1, hour=20, minute=3, second=0
    ).timestamp()
    exec_time = datetime(
        year=2021, month=1, day=1, hour=22, minute=13, second=0
    ).timestamp()
    test_result = test_subject.get_time_box_work(prev_exec_time, exec_time)
    assert test_result is not None, 'expect a result'
    assert len(test_result) == 2, 'wrong number of results'
    assert len(test_subject._json_cache._lookup) == 2, 'caching failed'
    test_entry = test_result.popleft()
    assert test_entry.entry_name == 'N20210101S0043.fits', 'wrong first file'
    assert test_entry.entry_ts == 1609535565.237183, 'wrong first timestamp'
    test_entry = test_result.popleft()
    assert test_entry.entry_name == 'N20210101S0042.fits', 'wrong 2nd file'
    assert test_entry.entry_ts == 1609535567.250666, 'wrong 2nd timestamp'

    # get nothing
    prev_exec_time = datetime(
        year=2019, month=1, day=1, hour=20, minute=3, second=0
    ).timestamp()
    exec_time = datetime(
        year=2019, month=2, day=1, hour=22, minute=13, second=0
    ).timestamp()
    test_result = test_subject.get_time_box_work(prev_exec_time, exec_time)
    assert test_result is not None, 'expect a result'
    assert len(test_result) == 0, 'wrong number of empty result list'
    assert len(test_subject._json_cache._lookup) == 0, 'flushing failed'
