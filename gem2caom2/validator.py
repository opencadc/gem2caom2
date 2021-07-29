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

import logging
import traceback

from datetime import datetime

from caom2pipe import manage_composable as mc
from gem2caom2 import external_metadata
from gem2caom2.util import COLLECTION, SCHEME

__all__ = ['GeminiValidator']


class GeminiValidator(mc.Validator):
    def __init__(self):
        super(GeminiValidator, self).__init__(
            source_name=COLLECTION,
            scheme=SCHEME,
            preview_suffix='_th.jpg',
        )
        self._gofr = external_metadata.get_gofr()
        self._max_date = datetime.fromtimestamp(
            self._gofr.get_max_timestamp()
        ).date()
        logging.error(f'max date is {self._max_date}')
        self._rejected = mc.Rejected(self._config.rejected_fqn)
        self._logger = logging.getLogger(__name__)

    def _date_file_name(self, file_name):
        """Extract the date from the file name. Use this value to compare
        against the maximum date at the source. This might avoid
        spurious reporting on files that are not yet listed at the source.
        It might not.
        """
        try:
            pattern = '%Y%m%d'
            splitter = 'S'
            index = 0
            candidate = file_name.split('_')[0]
            if candidate[0].isdigit():
                pattern = '%y%b%d'
                if candidate.count('.') == 2:
                    splitter = '.'
                else:
                    splitter = '_'
                    if candidate[2].isdigit():
                        pattern = '%Y%b%d'
            elif candidate[0].lower() == 'p':
                pattern = '%Y%b%d'
                index = 1
                splitter = '_'
            elif candidate.count('S') > 1:
                index = candidate.index('S') + 1
            elif candidate.count('N') == 1:
                if candidate.count('Q') == 1:
                    pattern = '%Y'
                    index = 2
                    # semesters
                    splitter = 'B'
                    if candidate.count('A') == 1:
                        splitter = 'A'
                else:
                    index = candidate.index('N') + 1
            elif candidate.count('TX') == 1:
                index = candidate.index('TX') + 2
            elif candidate[0] == 'r':
                pattern = '%y%b%d'
                index = 1
                splitter = '_'
            elif candidate[0] == 'c':
                pattern = '%Y%b%d'
                index = 1
                splitter = '_'
            elif candidate[0] == 'a' and candidate[1] == 'g':
                pattern = '%Y%b%d'
                index = 2
                splitter = '_'
            elif candidate[0] == 'G':
                pattern = '%Y'
                index = 2
                splitter = candidate[6]
            else:
                logging.warning(
                    f'Pretty sure this is an unexpected file name format '
                    f'{file_name}'
                )

            candidate = candidate[index:]
            if candidate.count('G') == 1:
                splitter = 'G'
            result = datetime.strptime(
                candidate.split(splitter)[0], pattern
            ).date()
        except ValueError as ex:
            self._logger.error(
                f'Do not understand date format in file name {file_name}'
            )
            result = None

        return result

    def _get_date_remove_set(self, coll, coll_name):
        remove = set()
        for f_name in coll:
            f_date = self._date_file_name(f_name)
            if f_date is not None and f_date > self._max_date:
                remove.add(f_name)
                logging.warning(f'Removing {f_name} from {coll_name}.')
        return remove

    def _filter_result(self):
        """Remove results from the lists returned from CADC data and metadata
        queries which are dated after the maximum timestamp in the input
        file, because there's no source information to which it's possible
        to compare those files.

        Not sure if the maximum timestamp of the supplied file list is a good
        termination time for comparisons, or not, but that's where the code
        is going to start.
        """
        remove = self._get_date_remove_set(self._source, 'source')
        self._source = list(set(self._source).difference(remove))
        remove = self._get_date_remove_set(
            self._destination_meta, 'destination meta'
        )
        self._destination_meta = list(
            set(self._destination_meta).difference(remove)
        )
        remove = self._get_date_remove_set(
            self._destination_data, 'destination data'
        )
        self._destination_data = list(
            set(self._destination_data).difference(remove)
        )

    def read_from_source(self):
        result = {}
        file_ids = self._gofr.name_list.keys()
        # begin with the assumption there's a preview file for every fits
        # file, but do not add entries for preview files that have been
        # identified as unavailable for some reason from archive.gemini.edu
        #
        # ignore fits files that have been rejected because of bad
        # metadata
        for file_id in file_ids:
            f_name = f'{file_id}.fits'
            if not self._rejected.is_bad_metadata(f_name):
                ts_s = self._gofr.get_timestamp(file_id)
                result[f_name] = ts_s

                f_name = f'{file_id}.jpg'
                if not self._rejected.is_no_preview(f_name):
                    result[f_name] = ts_s
        return result

    def write_todo(self):
        if len(self._source) == 0 and len(self._destination_data) == 0:
            logging.info(f'No entries to write to {self._config.work_fqn}')
        else:
            with open(self._config.work_fqn, 'w') as f:
                for entry in self._source:
                    f.write(f'{entry}\n')
                for entry in self._destination_data:
                    f.write(f'{entry}\n')


def validate():
    validator = GeminiValidator()
    validator.validate()
    validator.write_todo()


if __name__ == '__main__':
    import sys

    try:
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        validate()
        sys.exit(0)
    except Exception as e:
        logging.error(e)
        logging.debug(traceback.format_exc())
        sys.exit(-1)
