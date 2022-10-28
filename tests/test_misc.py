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

import shutil
import tempfile
import unittest
from scibacr import *
import pysam
import os


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

    def test_kmer_sliding(self):
        gen_kmer_sliding_window_ref(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta', 6,
                                    'data/kmer_6.h5', genes=['NC_003210.1:1000534-1001425', 'NC_003210.1:1008129-1008879'])

    def test_fasta_to_gff(self):
        convert_transcriptome_to_gff(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta',
                                     'data/gff.gff')

    def test_kmer_extraction(self):
        kmer = 'GAAGAT'
        kmer_h5 = h5py.File('data/kmer_6.h5', 'r')
        kmer_data = kmer_h5[kmer]
        print(kmer_data.keys())
        for k in kmer_data.keys():
            for position in kmer_data[k]:
                print(position)

    def test_dict_meta(self):
        read_dict = gen_mapping_gene_read_dict_meta(f'data/WT-K12_MG1655_genome.sorted.bam', f'data/ssu_all_r207_gr1400.gff')
        print(read_dict['RS_GCF_002942685.1~NZ_PTNW01000082.1'])
        assert '7bfc9ef7-f783-4ef0-a89a-fc56645687ff_Basecall_1D_template' in read_dict['RS_GCF_002942685.1~NZ_PTNW01000082.1']
        gen_training_h5py_position('data/WT-K12_MG1655_genome.sorted.bam', 'data/ssu_all_r207_gr1400.fna', read_dict,
                                   output_filename='data/training_meta.h5',
                                   min_coverage=20, max_coverage=100)

    def test_epinano(self):
        # python $EPINANO_HOME/Epinano_Variants.py -n 6 -R reference.fasta -b sample.reads.bam -s  $EPINANO_HOME/misc/sam2tsv.jar --type t
        os.system('python ')
    def test_counts_meta(self):
        df = count_reads_meta([f'data/WT-K12_MG1655_genome.sorted.bam'], f'data/ssu_all_r207_gr1400.gff', transcriptome=False)
        df.to_csv('data/metagenome.csv')

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
        ref = pysam.FastaFile(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta')
        pos = 'NC_003210.1:1546342-1546531'
        bam = pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
        write_msa_over_gene(pos, bam, ref, 'data/msa.fasta')

    def test_extract_quality(self):
        data = h5py.File(f'data/SRR13212638.h5', 'r')
        for k in data.keys(): #['NC_003210.1:1546342-1546531']:
            if k == 'NC_003210.1:1546342-1546531':
                for read in data[k]:
                    if 'SRR13212638.205:' in read:
                        print(k)
                        read_quality = data[k][read]
                        read_quality = [r for r in read_quality]
                        print(read_quality)
                        qual = [8,3,3,1,6,5,3,3,4,5,15,20,10,11,14,10,14,3,11,4,9,5,9,11,3,5,11,14,6,19,17,9,2,6,6,6,9,7,3,6,4,8,12,9,10,4,9,5,4,3,2,6,5,5,8,6,17,25,29,33,27,24,16,21,12,7,8,17,22,21,11,17,18,20,21,20,22,22,20,22,19,16,4,4,6,10,3,4,4,10,26,27.0,7.0,7.0,26.0,26.0,26.0,26.0,28.0,26.0,26.0,23.0,23.0,25.0,13.0,9.0,18.0,13.0,25.0,23.0,25.0,18.0,3.0,18.0,23.0,23.0,12.0,29.0,16.0,10.0,23.0,21.0,20.0,22.0,32.0,25.0,24.0,22.0,27.0,22.0,27.0,32.0,33.0,33.0,28.0,33.0,15.0,41.0,31.0,12.0,15.0,13.0,7.0,3.0,1.0,2.0,4.0,6.0,4.0,15.0,25.0,12.0,18.0,22.0,25.0,25.0,29.0,33.0,38.0,32.0,4.0,10.0,29.0,12.0,35.0,8.0,17.0,13.0,2.0,3.0,18.0,21.0,18.0,14.0,17.0,16.0,8.0,12.0,25.0,34.0,21.0,22.0,24.0,21.0,10.0,11.0,7.0,9.0,17.0,24.0,24.0,13.0,15.0,18.0,30.0,39.0,38.0,31.0,42.0,37.0,28.0,27.0,29.0,30.0,31.0,33.0,36.0,29.0,38.0,36.0,36.0,37.0,33.0,33.0,34.0,32.0,41.0,34.0,38.0,27.0,27.0,29.0,38.0,41.0,37.0,24.0,24.0,25.0,20.0,26.0,15.0,21.0,20.0,16.0]
                        print(qual)
                        assert np.sum(read_quality) == np.sum(qual, dtype=int)
                        break
            break

    def test_training_from_h5py(self):
        # Check that we can get the same as we got in the csv
        gene = 'NC_003210.1:1546342-1546531'  #
        read = 'SRR13212638.205'
        data = h5py.File(f'data/training.h5', 'r')
        print(data[gene][read]['qual'][0]) # 8
        seq = [chr(s) for s in data[gene][read]['seq']]
        print(''.join([chr(s)for s in data[gene][read]['seq']]))
        print(data[gene][read]['seq'])
        qual = [s for s in data[gene][read]['qual']]
        print(len(seq), len(qual))

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

    def test_count_reads(self):
        read_df = count_reads([f'data/SRR13212638.sorted.bam'], f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genes-RNA.gff')
        read_df.to_csv(f'data/counts.csv', index=False)

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
                        seq, qual, ref, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str, read.query_qualities)
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
        #mm = get_mismatches(samfile, ref, label, 1546342, 1546531)  # Don't think this works correctly awks
        #print(mm)

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
        # #Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,de
        #  python Epinano_Predict.py  --predict /Users/ariane/Documents/code/scibacr/tests/data/SRR13212638.sorted.plus_strand.per.site.csv --out_prefix some_sample.modification --columns=7,8,10 --model models/rrach.q3.mis3.del3.linear.dump
        #
        # /Users/ariane/Documents/code/scibacr/tests/data/SRR13212638.sorted.bam
        # /Users/ariane/Documents/code/scibacr/tests/data/SRR13212638.sorted.plus_strand.per.site.csv
        # /Users/ariane/Documents/code/scibacr/tests/data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta
        samfile = pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
        ref = pysam.FastaFile(f'data/SRR13212638_GCF_000196035.1_ASM19603v1_genomic_transcripts.fasta')
        pos = 'NC_003210.1:1013474-1014527'
        reads = []
        for read in samfile.fetch(pos):
            if read.is_reverse == False:
                reads.append(read)
        read = reads[0]
        ref_str = ref[pos]
        seq, ref, ins, qual = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str, read.query_qualities)
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

