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
import pandas as pd
import numpy as np
import pysam
from tqdm import tqdm


def dedup_gff(gff_file: str, output_file=None, feature_filter=None, source_filter=None) -> pd.DataFrame:
    """
    Removes duplicate entries from a gff file. This can be needed when mapping reads to transcriptome.

    Fields: https://www.ensembl.org/info/website/upload/gff.html?redirect=no
    Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position* of the feature, with sequence numbering starting at 1.
    end - End position* of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.


    Parameters
    ----------
    gff_file: gff file
    output_file: optional, the output path
    feature_filter: optional, list of gene/RNA names to filter, e.g. ['gene', 'tRNA', 'ncRNA', 'rRNA']
    source_filter: optional, list of sources you want to filter e.g. ['RefSeq']

    Returns
    -------

    """
    rows = []
    output_file = output_file or gff_file.replace('.gff', '_simple.gff')
    with open(f'{gff_file}', 'r') as gff:
        for line in gff:
            if line[0] == '#':
                continue
            else:
                line = line.split('\t')
                line[-1] = line[-1].strip()
                rows.append([f'{line[0]}:{line[3]}-{line[4]}'] + line)
    df = pd.DataFrame(rows)
    if source_filter is not None:
        df = df[df[2].isin(source_filter)]
    # OK so let's just keep the genes and the tRNA, rRNA, and ncRNAs
    if feature_filter is not None:
        df = df[df[3].isin(feature_filter)]
    # Sort by '3' and drop duplicates (sort by feature)
    df = df.sort_values(by=[0, 3])
    df = df.drop_duplicates(subset=[0], keep='first')
    # Save to tsv without the header and also drop that first column (i.e. save rows 1 - 9)
    df = df.drop(columns=[0])
    df.to_csv(f'{output_file}', sep='\t', index=False, header=False)
    return df


def gen_training_data(bam, ref, gff, min_coverage=20, max_coverage=100):
    """

    Parameters
    ----------
    bam
    ref
    gff
    min_coverage
    max_coverage

    Returns
    -------

    """
    return


import math
def chunk_data(x, y, chunk_size=7): # Use 7 as this means we can relate to eligos2
    # https://stackoverflow.com/questions/13673060/split-string-into-strings-by-length
    chunks, chunk_size = len(x), chunk_size
    return [x[i:i+chunk_size] for i in range(0, chunks, chunk_size)], [y[i:i+chunk_size] for i in range(0, chunks, chunk_size)]


def gen_training(bam, ref, gff, output_filename=None, kmer_len=7, min_coverage=20, max_coverage=1000):
    bam = pysam.AlignmentFile(bam, "rb")
    fasta = pysam.FastaFile(ref)
    gff = pd.read_csv(gff, header=None, sep='\t')

    positions = []
    pos_to_info = {}
    for line in gff.values:
        positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
        # Also keep track of the gene information
        pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()

    encodings = []
    for pos in tqdm(positions):
        reads = []
        for read in bam.fetch(pos):
            reads.append(read)
        read_num = 0
        if len(reads) > min_coverage:
            ref_str = fasta[pos]

            for ri, read in enumerate(reads):
                if read.query_sequence is not None:
                    seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                               read.query_qualities)
                    name = read.query_name
                    # Make it totally align
                    seq = "*" * read.reference_start + seq + "-" * (len(ref_str) - (read.reference_start + len(seq)))
                    qual = ([None] * read.reference_start) + qual + (
                                [None] * (len(ref_str) - (read.reference_start + len(seq))))
                    # Chunk the data...
                    chunked_seq, chunked_qual = chunk_data(seq, qual, kmer_len)

                    # Now we want to one hot encode the seq but add in the qual info...
                    for i, s in enumerate(chunked_seq):
                        enc = one_hot_gencode(s, chunked_qual[i])
                        encodings.append(
                            [name, pos, i * kmer_len] + list(
                                enc.flatten()))  # Want to flatten it so that we can get it easy
                    read_num += 1
                    if read_num > max_coverage:
                        break
    bam.close()
    df = pd.DataFrame(encodings)
    df = df.dropna()
    if output_filename is not None:
        df.to_csv(output_filename, index=False)
    return df


def gen_coverage_over_genes(bam, ref, gff, min_coverage=20, max_coverage=1000):
    """

    Parameters
    ----------
    bam
    ref
    gff
    min_coverage
    max_coverage

    Returns
    -------

    """
    bam = pysam.AlignmentFile(bam, "rb")
    fasta = pysam.FastaFile(ref)
    gff = pd.read_csv(gff, header=None, sep='\t')

    positions = []
    pos_to_info = {}
    for line in gff.values:
        positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
        # Also keep track of the gene information
        pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()

    gene_rows = []
    for pos in tqdm(positions):
        reads = []
        for read in bam.fetch(pos):
            reads.append(read)
        read_num = 0
        if len(reads) > min_coverage:
            ref_str = fasta[pos]
            rows = []
            for ri, read in enumerate(reads):
                if read.query_sequence is not None:
                    seq, ref, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str)
                    row = [read.query_name]
                    # Make it totally align
                    seq = "*" * read.reference_start + seq + "-" * (
                            len(ref_str) - (read.reference_start + len(seq)))
                    row += list(np.array(list(seq)))
                    rows.append(row)
                    read_num += 1
                    if read_num > max_coverage:
                        break
            rows = np.array(rows)
            try:
                sumarised = []  # np.zeros(len(ref_str)*4)
                for i in range(0, len(ref_str)):
                    sumarised.append(len(rows[rows[:, i] == 'A']))
                    sumarised.append(len(rows[rows[:, i] == 'T']))
                    sumarised.append(len(rows[rows[:, i] == 'G']))
                    sumarised.append(len(rows[rows[:, i] == 'C']))
                    sumarised.append(len(rows[rows[:, i] == '-']))  # Missing content (i.e. a deletion)
                sumarised_arr = np.array(sumarised)
                gene_rows.append([pos, pos_to_info[pos], np.mean(sumarised_arr), np.max(sumarised_arr),
                                  np.var(sumarised_arr)] + sumarised)
            except:
                x = 1
    df = pd.DataFrame(gene_rows)  # Can't give it columns since there will be differing column lengths
    bam.close()
    return df


def write_msa_over_gene(gene_location: str, bam, ref, output_filename=None):
    """
    Makes a pileup over a gene and saves it as a MSA so one can view it in AliView.

    Rows are the reads, columns are the columns in the reference. Insertions are ignored.
    Parameters
    ----------
    gene_location: location (chr:start-end) note start might need to be -1 depending on how the transcriptome was created
    bam: bam file read in by pysam, pysam.AlignmentFile(f'data/SRR13212638.sorted.bam', "rb")
    ref: fasta file read in by pysam, pysam.FastaFile(sam/bam file)
    output_filename

    Returns
    -------

    """
    reads = []
    for read in bam.fetch(gene_location):
        reads.append(read)
    output_filename = output_filename or 'reads_msa_{gene_location}.fasta'
    with open(output_filename, 'w') as fout:
        ref_str = ref[gene_location]
        fout.write(f'>Reference\n')
        fout.write(ref_str + '\n')
        for read in reads:
            seq, ref, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str)
            # Make it totally align
            seq = "-" * read.reference_start + seq + "-" * (len(ref_str) - (read.reference_start + len(seq)))
            fout.write(f'>{read.query_name}\n')
            fout.write(f'{seq}\n')


def one_hot_gencode(seq, qual):
    """
    One hot encode but put in the quality info of that position
    Parameters
    ----------
    seq
    qual

    Returns
    -------

    """
    encoded = np.zeros((4, len(seq)))
    for i, n in enumerate(seq):
        if n == 'A':
            encoded[0, i] = qual[i]  # We use qual to get the value
        elif n == 'T':
            encoded[1, i] = qual[i]  # We use qual to get the value
        elif n == 'G':
            encoded[2, i] = qual[i]  # We use qual to get the value
        elif n == 'C':
            encoded[3, i] = qual[i]  # We use qual to get the value
        else:
            encoded[0, i] = None  # This will be dropped later ToDo: think of a better way to incorperate deletions...
    return encoded


def alignment_from_cigar(cigar: str, alignment: str, ref: str, query_qualities: list) -> tuple[str, str, list, list]:
    """
    Generate the alignment from the cigar string.
    Operation	Description	Consumes query	Consumes reference
    0 M	alignment match (can be a sequence match or mismatch)	yes	yes
    1 I	insertion to the reference	yes	no
    2 D	deletion from the reference	no	yes
    3 N	skipped region from the reference	no	yes
    4 S	soft clipping (clipped sequences present in SEQ)	yes	no
    5 H	hard clipping (clipped sequences NOT present in SEQ)	no	no
    6 P	padding (silent deletion from padded reference)	no	no
    7 =	sequence match	yes	yes
    8 X	sequence mismatch	yes	yes

    Parameters
    ----------
    cigar: 49S9M1I24M2D9M2I23M1I3M1D13M1D10M1D13M3D12M1D23M1
    alignment: GCTGATCACAACGAGAGCTCTCGTTGCTCATTACCCCTAAGGAACTCAAATGACGGTTAAAAACTTGTTTTGCT
    ref: reference string (as above but from reference)

    Returns
    -------

    """
    new_seq = ''
    ref_seq = ''
    qual = []
    inserts = []
    pos = 0
    ref_pos = 0
    for op, op_len in cigar:
        if op == 0:  # alignment match (can be a sequence match or mismatch)
            new_seq += alignment[pos:pos+op_len]
            qual += query_qualities[pos:pos + op_len]

            ref_seq += ref[ref_pos:ref_pos + op_len]
            pos += op_len
            ref_pos += op_len
        elif op == 1:  # insertion to the reference
            inserts.append(alignment[pos - 1:pos+op_len])
            pos += op_len
        elif op == 2:  # deletion from the reference
            new_seq += '-'*op_len
            qual += [-1]*op_len
            ref_seq += ref[ref_pos:ref_pos + op_len]
            ref_pos += op_len
        elif op == 3:  # skipped region from the reference
            new_seq += '*'*op_len
            qual += [-2] * op_len
            ref_pos += op_len
        elif op == 4:  # soft clipping (clipped sequences present in SEQ)
            inserts.append(alignment[pos:pos+op_len])
            pos += op_len
        elif op == 5:  # hard clipping (clipped sequences NOT present in SEQ)
            continue
        elif op == 6:  # padding (silent deletion from padded reference)
            continue
        elif op == 7:  # sequence mismatch
            new_seq += alignment[pos:pos + op_len]
            ref_seq += ref[ref_pos:ref_pos + op_len]
            qual += query_qualities[pos:pos + op_len]
            pos += op_len
            ref_pos += op_len
    return new_seq, ref_seq, qual, inserts


def alignment_from_cigar_inc_inserts(cigar: str, alignment: str, ref: str) -> tuple[str, str, list]:
    """
    Generate the alignment from the cigar string.
    Operation	Description	Consumes query	Consumes reference
    0 M	alignment match (can be a sequence match or mismatch)	yes	yes
    1 I	insertion to the reference	yes	no
    2 D	deletion from the reference	no	yes
    3 N	skipped region from the reference	no	yes
    4 S	soft clipping (clipped sequences present in SEQ)	yes	no
    5 H	hard clipping (clipped sequences NOT present in SEQ)	no	no
    6 P	padding (silent deletion from padded reference)	no	no
    7 =	sequence match	yes	yes
    8 X	sequence mismatch	yes	yes

    Parameters
    ----------
    cigar: 49S9M1I24M2D9M2I23M1I3M1D13M1D10M1D13M3D12M1D23M1
    alignment: GCTGATCACAACGAGAGCTCTCGTTGCTCATTACCCCTAAGGAACTCAAATGACGGTTAAAAACTTGTTTTGCT
    ref: reference string (as above but from reference)

    Returns
    -------

    """
    new_seq = ''
    ref_seq = ''
    pos = 0
    ref_pos = 0
    for op, op_len in cigar:
        if op == 0:  # alignment match (can be a sequence match or mismatch)
            new_seq += alignment[pos:pos+op_len]
            ref_seq += ref[ref_pos:ref_pos + op_len]
            pos += op_len
            ref_pos += op_len
        elif op == 1:  # insertion to the reference
            new_seq += alignment[pos:pos + op_len]
            ref_seq += '^'*op_len
            pos += op_len
        elif op == 2:  # deletion from the reference
            new_seq += '-'*op_len
            ref_seq += ref[ref_pos:ref_pos + op_len]
            ref_pos += op_len
        elif op == 3:  # skipped region from the reference
            new_seq += '*'*op_len
            ref_pos += op_len
        elif op == 4:  # soft clipping (clipped sequences present in SEQ)
            new_seq += alignment[pos:pos + op_len]
            ref_seq += '>' * op_len
            pos += op_len
        elif op == 5:  # hard clipping (clipped sequences NOT present in SEQ)
            continue
        elif op == 6:  # padding (silent deletion from padded reference)
            continue
        elif op == 7:  # sequence mismatch
            new_seq += alignment[pos:pos + op_len]
            ref_seq += ref[ref_pos:ref_pos + op_len]
            pos += op_len
            ref_pos += op_len
    return new_seq, ref_seq