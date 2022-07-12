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
import pysam
import shutil
import tempfile
import unittest
import numpy as np
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

    def test_reverse_alignment(self):
        """
        Seems to work the same? Looks like they reverse it and then just keep the tag for future references.
        """
        samfile = pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
        ref = pysam.FastaFile(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta')
        pos = 'NC_003210.1:1029207-1029819'
        reads = []
        for read in samfile.fetch(pos):
            if read.is_reverse == True:
                reads.append(read)
        read = reads[0]
        ref_str = ref[pos]
        seq, ref, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str)
        # Check the similarity
        sim = 0
        for i, s in enumerate(seq):
            if s == ref[i] and s != '-':
                sim += 1
        print(sim, len(seq))
        print("new sequence & aligned reference")
        print(seq)
        print(ref)
        # We see that this is exactly the same :D
        print(ref_str[read.reference_start:read.reference_start + len(read.query_sequence)])

        print(read.get_aligned_pairs())
        print("Old query_alignment_sequence & old reference")
        print(read.query_alignment_sequence)
        ref_pos = read.get_reference_positions()
        ref_arr = np.array(list(ref_str))
        ref_str_from_aln = ''.join(ref_arr[ref_pos])
        print(ref_str_from_aln)

    def test_alignment(self):
        """
        SRR13212638.113019	0	NC_003210.1:1013474-1014527	1	60	19S6M1D3M1D9M1I14M1D15M3D28M1D12M1D15M1D35M1I22M1D41M2D20M1D4M1D43M1I15M2D20M1I6M1I20M1I22M2D7M2D13M1I7M1D48M1I19M1D10M1D3M1D10M1D2M1D11M1S	*00	ACACTAAAGGAGGAATTTTTTGTCAACGATCAAATACAAAAAGAAACAATGCAAAAATTAAAGAATTATCTATCCCAAGCCCAACTGGCAATACAGGAAAATCATCGAAAACTGGCAAAAGATTGGATGCTTCAAAAATCCCTTACCATTTAAACAACAAAAGGTGGATTAATCGCAACGTTCCTGGGAAAGACGACACGAAACACCGTATGCTAACCGCCTGTGGATACGCTGGGTGCAAGCTTGAGAAATTAAAGCAGACGGACGATTACTTTTAACCTTAATTGGAGGGCTACCGTTTTAACGATTGAAGGCGAATATTGTCACCATTGGAGACAAGTGATGGAGATTTCGTATAGTGGCACCATTTTAATGTCAAACTTGTCCATGTTTATAAAAGATGTTGGACTGCTGAGCGTAATGACAAAAATATGGAAGTTCGCTAAGATGTAAAAAGCATTAGATGCAGACAAGTTCGTGCGTTGAATAGAAGTGGGATTTTGTTTCG	%'$'0010/C=58DB>A;>;4<6&'2461&/,87<>><**+=>9==:2005%/3=>>D:;-*+>94?F;2>=:-;BE?+6;:/../5<325084134<78<>;>65AFGD=A>;6:5:9FAC442<BB@,(7.1:<B@5&:.>4:<?:-*-62$-4>=?=>?B?D=>7.26//%2,3@:?%/-.%%7<5-3=:04/.+35,,)+$%#2,)-82.+'(,2.##(*3?110*+2++++%$""#/##%$,3&68999<<B11?)).>01G**))))&%%('($().'/078'53:5H),1=FKCA<0.3?:<:QF8<:@@4(-607+$$.0*/3.*)7>;69D@BDE?=FAB>=9->7<27<DHE:6;681?@?;28>57::93):?3-64372*560=?;=<0+(*::9@7?<58D<=D=+553B:3:9:99:502C559887735$/<25<114435<4&'0$<81<:345357700064<):%*5,,*(>&'7&*2>&&&&&8$$$$6	NM:i:69	ms:i:288	AS:i:285	nn:i:0	ts:A:+	tp:A:P	cm:i:29	s1:i:202	s2:i:0	de:f:0.124rl:i:0
        """
        samfile = pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
        ref = pysam.FastaFile(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta')
        pos = 'NC_003210.1:1013474-1014527'
        reads = []
        for read in samfile.fetch(pos):
            if read.is_reverse == False:
                reads.append(read)
        read = reads[0]
        ref_str = ref[pos]
        seq, ref, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str)
        # Check the similarity
        sim = 0
        for i, s in enumerate(seq):
            if s == ref[i] and s != '-':
                sim += 1
        print(sim, len(seq))
        # Similarity should be the number of matches
        match_num = 0
        mismatch = 0
        dels = 0
        for op, oplen in read.cigartuples:
            if op == 0:  # alignment match (can be a sequence match or mismatch)
                match_num += oplen
            elif op == 7:
                mismatch += oplen
            elif op == 2:
                dels += oplen
        print(match_num, mismatch, dels)
        assert match_num + mismatch + dels == len(seq)
        assert len([c for c in seq if c == '-']) == dels
        print("new sequence & aligned reference")
        print(seq)
        print(ref)
        # We see that this is exactly the same :D
        print(ref_str[read.reference_start:read.reference_start + len(read.query_sequence)])

        print(read.get_aligned_pairs())
        print("Old query_alignment_sequence & old reference")
        print(read.query_alignment_sequence)
        ref_pos = read.get_reference_positions()
        ref_arr = np.array(list(ref_str))
        ref_str_from_aln = ''.join(ref_arr[ref_pos])
        print(ref_str_from_aln)
        old_sim = 0
        for i, s in enumerate(read.query_alignment_sequence):
            if i < len(ref_str_from_aln):
                if s == ref_str_from_aln[i] and s != '-':
                    old_sim += 1
        # We see that there is a very evident mis-match...they don't add in the deletions and include the insertions
        # which means we can't identify the differences easily for every position in the reference...
        print(old_sim, len(read.query_alignment_sequence), len(ref_str_from_aln))

        print("Input query_sequence & insertions")
        print(read.query_sequence)
        print(ins)

        # Looks correct :)
        print("As above but including the inserted sequences from the read")
        s = alignment_from_cigar_inc_inserts(read.cigartuples, read.query_sequence, ref_str)
        print(s[0])
        # It removes the clipping so add this in to look at it
        print('GCCACCCCGAAAAGGAGCGTTTTGTGCTATTT' + read.query_alignment_sequence)
        print(s[1])
        print(read.query_sequence)
        #print(read.get_reference_sequence())
        # reads.is_reverse
        # reads.query_name
        # reads.seq
        # reads.mapping_quality
        # reads.get_cigar_stats()
        # reads.query_qualities  # Sequence quality
        # reads.query_sequence
        # reads.query_alignment_end
        # reads.query_alignment_length
        # reads.query_alignment_qualities
        # reads.query_alignment_sequence

