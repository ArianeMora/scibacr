###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import pandas as pd
import shutil
import tempfile
import unittest

from scibacr import *

class TestExample(unittest.TestCase):

    def setUp(self):
        # Flag to set data to be local so we don't have to download them repeatedly. ToDo: Remove when publishing.
        self.local = True

        if self.local:
            THIS_DIR = os.path.dirname(os.path.abspath(__file__))
            self.tmp_dir = os.path.join(THIS_DIR, 'data/tmp/')
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='EXAMPLE_PROJECT_tmp_')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_dedup_gff(self):
        gff = 'data/GCF_000026545.1_ASM2654v1_genomic.gff'
        old_gff = pd.read_csv(gff, sep='\t', header=None)
        new_gff = dedup_gff(gff, None, None, None)
        # Check it actually filtered
        assert len(new_gff) < len(old_gff)
        # Check they have the same IDs (i.e. chr, start, end)
        old_ids = []
        for line in old_gff.values:
            try:
                old_ids.append(f'{line[0]}:{int(line[3])}-{int(line[4])}')
            except:
                print(line)
        new_ids = []
        for line in new_gff.values:
            try:
                new_ids.append(f'{line[0]}:{int(line[3])}-{int(line[4])}')
            except:
                print(line)
        # i.e. check everything the same when we have no filters!
        assert len(set(old_ids) & set(new_ids)) == len(set(old_ids))
        # Check filters
        new_gff = dedup_gff(gff, None, None, ['RefSeq'])
        assert len(new_gff) < len(old_gff)
        # Check only refseq in the column
        f = list(set(new_gff.values[:, 1]))
        assert f == ['RefSeq']

        # Check filters
        new_gff = dedup_gff(gff, None, ['tRNA'], None)
        assert len(new_gff) < len(old_gff)
        # Check only refseq in the column
        f = list(set(new_gff.values[:, 2]))
        assert f == ['tRNA']



