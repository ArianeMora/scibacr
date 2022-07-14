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
import pysam
import re, os, sys


def filter_reads(reads, min_mapq=0):
    #     reads = [read for read in reads if not read.alignment.is_supplementary and not read.alignment.is_unmapped and not read.alignment.is_duplicate and not read.is_refskip]
    #     if min_mapq > 0:
    #         reads = [read for read in reads if read.alignment.mapq >= min_mapq]
    return reads


def divided(x, y):
    if x == 0 or y == 0:
        return float(0)
    else:
        return float(x) / float(y)


##### Decision modules
def rvStrand(reverse, is_reverse, bam_type=None):
    if bam_type == 'cdna':
        return True
    elif reverse == True and is_reverse:
        return True
    elif reverse == False and not is_reverse:
        return True
    else:
        return False


bam_type = None

def get_mismatches(bamData, faidx, chrom, start, end, name=None, bam_type=None, reverse=False, max_depth=10000):
    homopolymer = re.compile(r'(A{4,}|T{4,}|C{4,}|G{4,})', re.I)
    outTmp = []

    for pileupcolumn in bamData.pileup(chrom):
        read_counts = 0
        reads = filter_reads(pileupcolumn.pileups)
        name = chrom
        # reads_nodel = [read for read in reads if not read.is_del]
        start_b = pileupcolumn.pos
        end_b = start_b + 1
        start_loc = start + start_b
        end_loc = start_loc + 1

        ref = faidx[chrom][start_b:end_b].upper()
        kmer5 = faidx[chrom][start_b - 2:end_b + 2].upper()
        kmer7 = faidx[chrom][start_b - 2:end_b + 4].upper()
        strand = '-' if reverse else '+'

        reads_o = len([read for read in reads if rvStrand(reverse, read.alignment.is_reverse, bam_type)])
        matches_o, substitutions_o, substitutions_wo_ins_o = (0, 0, 0)
        basesCounts = {"A": 0, "T": 0, "C": 0, "G": 0}

        cdnaMajorAllele = False
        for read in reads:
            read_counts += 1
            ## count matches, substitutions
            if not read.is_del:
                if rvStrand(reverse, read.alignment.is_reverse, bam_type):
                    basesCounts[read.alignment.seq[read.query_position]] += 1
                    if cdnaMajorAllele:
                        if read.alignment.seq[read.query_position] == cdnaMajorAllele:
                            matches_o += 1
                        else:
                            substitutions_o += 1
                            if not read.indel > 0:
                                substitutions_wo_ins_o += 1
                    else:
                        if read.alignment.seq[read.query_position] == ref:
                            matches_o += 1
                        else:
                            substitutions_o += 1
                            if not read.indel > 0:
                                substitutions_wo_ins_o += 1

        ## Major allel
        majAllel = max(basesCounts, key=lambda k: basesCounts[k])
        majAllelFreq = divided(basesCounts[majAllel], float(sum(basesCounts.values())))
        ## homopolymer count
        homoseq = "--"
        test_start_home = start_b - 2 if start_b - 2 > 0 else 0
        homo = homopolymer.search(faidx[chrom][test_start_home:start_b + 4])
        if (homo):
            homoseq = homo.group(1).upper()

        deletions_o = len([read for read in reads
                           if read.is_del and rvStrand(reverse, read.alignment.is_reverse, bam_type)])
        insertions_o = len([read for read in reads
                            if read.indel > 0 and rvStrand(reverse, read.alignment.is_reverse, bam_type)])
        error_o = substitutions_wo_ins_o + insertions_o + deletions_o
        o_ref = ''#complement(ref) if reverse else ref
        o_kmer5 = ''#reversecomplement(kmer5) if reverse else kmer5
        o_homoseq = ''#reversecomplement(homoseq) if reverse else homoseq
        o_kmer7 = ''#reversecomplement(kmer7) if reverse else kmer7
        o_majAllel = ''#reversecomplement(majAllel) if reverse else majAllel
        ## for model matching
        m_kmer5 = kmer5
        m_kmer7 = kmer7
        name = None
        outTmp.append([chrom, start_loc, end_loc, strand, name, o_ref,
                       reads_o, matches_o, error_o,
                       substitutions_o, deletions_o, insertions_o,
                       divided(error_o, reads_o),
                       divided(substitutions_o, reads_o),
                       divided(deletions_o, reads_o),
                       divided(insertions_o, reads_o),
                       o_homoseq, o_kmer5, o_majAllel, majAllelFreq, o_kmer7,
                       m_kmer5, m_kmer7
                       ])
    return outTmp


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

    def test_msa(self):
        samfile = pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
        ref = pysam.FastaFile(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta')
        pos = 'NC_003210.1:1546342-1546531'
        write_msa_over_gene(pos, samfile, ref, 'data/msa.fasta')

    def test_cov_df(self):
        df = gen_coverage_over_genes(f'data/SRR13212638.sorted.bam',
                                     f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta',
                                     'data/SRR13212638_GCF_000196035.1_ASM19603v1_genes-RNA.gff', min_coverage=0,
                                     max_coverage=20)
        print(df.head())
        df.to_csv('data/coverage.csv')

    def test_training(self):
        df = gen_training(f'data/SRR13212638.sorted.bam',
                          f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta',
                          f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genes-RNA.gff')
        print(df.head())

    def test_training_gen(self):
        samfile = pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
        ref = pysam.FastaFile(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta')
        gff = pd.read_csv(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genes-RNA.gff', header=None, sep='\t')
        reads = []
        pos = 'NC_003210.1:1546342-1546531'
        name = 'SRR13212638'
        kmer_len = 7
        for read in samfile.fetch(pos):
            reads.append(read)
        if len(reads) > 100:
            ref_str = ref[pos]
            encodings = []
            for read in reads:
                seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                           read.query_qualities)
                row = [read.query_name]
                # Make it totally align
                seq = "*" * read.reference_start + seq + "-" * (len(ref_str) - (read.reference_start + len(seq)))
                qual = ([None] * read.reference_start) + qual + ([None] * (len(ref_str) - (read.reference_start + len(seq))))
                row += list(np.array(list(seq)))
                # Chunk the data...
                chunked_seq, chunked_qual = chunk_data(seq, qual, kmer_len)

                # Now we want to one hot encode the seq but add in the qual info...
                for i, s in enumerate(chunked_seq):
                    enc = one_hot_gencode(s, chunked_qual[i])
                    encodings.append([name, pos, i*kmer_len] + list(enc.flatten()))  # Want to flatten it so that we can get it easy
            df = pd.DataFrame(encodings)
            df = df.dropna()
            df.to_csv('data/training.csv', index=False)

            # Check it's the same as the pileup...

    def test_cov(self):
        samfile = pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
        ref = pysam.FastaFile(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta')
        gff = pd.read_csv(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genes-RNA.gff', header=None, sep='\t')
        reads = []
        pos = 'NC_003210.1:1546342-1546531'
        for read in samfile.fetch(pos):
            reads.append(read)
        if len(reads) > 100:
            ref_str = ref[pos]
            read_info = []
            rows = []
            for read in reads:
                seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                           read.query_qualities)
                row = [read.query_name]
                # Make it totally align
                seq = "*" * read.reference_start + seq + "-" * (
                            len(ref_str) - (read.reference_start + len(seq)))
                row += list(np.array(list(seq)))
                read_info.append({
                    'start': read.reference_start,
                    'end': read.reference_end,
                    'seq': seq,
                    'ins': ins,
                    'name': read.query_name
                })
                rows.append(row)
            rows = np.array(rows)
            sumarised = [] #np.zeros(len(ref_str)*4)
            for i in range(0, len(ref_str)):
                sumarised.append(len(rows[rows[:, i] == 'A']))
                sumarised.append(len(rows[rows[:, i] == 'T']))
                sumarised.append(len(rows[rows[:, i] == 'G']))
                sumarised.append(len(rows[rows[:, i] == 'C']))
                sumarised.append(len(rows[rows[:, i] == '-']))  # Missing content (i.e. a deletion)
            df = pd.DataFrame(rows, columns=["Name"] + [f'{c}-{i}' for i, c in enumerate(ref_str)] )#list(np.array(list(ref_str))))
            print(df)
            # Get value counts for each columns
            # Now we can just summarise that
            print(df['C-0'].value_counts())
            print(df['A-3'].value_counts())
            sumarised = np.array(sumarised)
            print(sumarised.reshape((len(ref_str), 5)))
            sumarised = sumarised.reshape((len(ref_str), 5))
            assert sumarised[1][3] == 77  # Same as in the value counts
            assert sumarised[4][2] == 3  # Same as in the value counts
            assert sumarised[4][0] == 83  # Should be
            assert sumarised[4][0] + sumarised[4][2] == 86  #
            # Check it's the same as the pileup...
            """
            NC_003210.1:1546342-1546531     1       C       77      ^],^],^\,^],^],^\,^],^I,^],^],^],^],^Z,^=,^],^],^2,^],^H,^H,^],^M,^L,^],^],^],^L,^],^],^],^],^Z,^],^3,^.,^A,^],^Y,^],^],^],^],^],^],^],^],^],^<,^X,^],^],^],^C,^],^],^],^],^8,^],^],^$,^],^=,^H,^N,^],^],^],^*,^],^],^],^],^],^@,^],^], !!!!!!!!!!!=!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!>!!!!!!!!!<!!!!!!!!!!!!!!!!!!@!!<!
            NC_003210.1:1546342-1546531     2       C       80      ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,^%,^],^3,  !%,!!!!,)!(4!,$!!!!!$$$"!'!!$!!!!-!$!!!!!!!<!!!!(!!!!C!$!!!,!!!!!!!!!!!!:!!1(!!(
            NC_003210.1:1546342-1546531     3       A       81      ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,^],     !71!!/!00!*0!,+!!!!!5(/#!,!!1!!!!5!6!!!!!!!9!!&!0!!!!3!:!!!4%!!!!!!!'!!!2!$*+&2+!
            NC_003210.1:1546342-1546531     4       A       86      ,,,,,,g,,,,,,,,,,,,g,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,-1t,,,g,,,,,,,,,,,,,,^;,^I,^W,^],^],     !-./!%!-5!'/!1+!!!!!0,&#!2!!1!!!!4!)!!!!!!!7!!+!5!!!!2!:!!!5(!!!!!!!*!!!*!+..)0$!!(!,%
            NC_003210.1:1546342-1546531     5       T       87      ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,,^2,       !$03!.!/7!#0!2+!!!!!5%$#!9!!1!!!!4!9!!!!!!!8
            """

    def test_pileup(self):
        """
        For a given gene we want to get the letters that map to each position....
        """
        samfile = pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
        ref = pysam.FastaFile(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta')
        gff = pd.read_csv(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genes-RNA.gff', header=None, sep='\t')
        positions = []
        for line in gff.values:
            try:
                positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
            except:
                print(line)
        for pos in positions:
            reads = []
            for read in samfile.fetch(pos):
                reads.append(read)
            if len(reads) > 100:
                with open(f'data/reads_msa_{pos}.fasta', 'w') as fout:
                    ref_str = ref[pos]
                    read_info = []
                    rows = []
                    fout.write(f'>Reference\n')
                    fout.write(ref_str + '\n')
                    for read in reads:
                        seq, ref, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str, read.query_qualities)
                        row = [read.query_name]
                        # Make it totally align
                        seq = "-"*read.reference_start + seq + "-"*(len(ref_str) - (read.reference_start + len(seq)))
                        fout.write(f'>{read.query_name}\n')
                        fout.write(f'{seq}\n')
                        row += list(np.array(list(seq)))
                        read_info.append({
                            'start': read.reference_start,
                            'end': read.reference_end,
                            'seq': seq,
                            'ins': ins,
                            'name': read.query_name
                        })
                        rows.append(row)
                    df = pd.DataFrame(rows, columns=["Name"] + list(np.array(list(ref_str))))
                    print(df)
                    break
            # except:
            #     print(pos)
        # NC_003210.1:1546342-1546531
        # Check pileup column implementation
        label = 'NC_003210.1:1546342-1546531'
        ref = pysam.FastaFile(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta')
        mm = get_mismatches(samfile, ref, label, 1546342, 1546531)  # Don't think this works correctly awks
        print(mm)

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
        seq, ref, ins, qual = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                   read.query_qualities)
        print("XXXXX", len(read.query_qualities), len(read.query_sequence))
        # Check the similarity
        sim = 0
        for i, s in enumerate(seq):
            if s == ref[i] and s != '-':
                sim += 1
        print(sim, len(seq))
        print("new sequence & aligned reference")
        print(seq)
        print(ref)
        print(''.join([str(i) for i in ins]))
        # We see that this is exactly the same :D
        print(ref_str[read.reference_start:read.reference_end])

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
        seq, ref, ins, qual = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str, read.query_alignment_qualities)
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
        s = alignment_from_cigar_inc_inserts(read.cigartuples, read.query_sequence, ref_str, read.query_alignment_qualities)
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

