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
from collections import defaultdict
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import numpy as np
import random
import pysam
import h5py
from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from wordcloud import WordCloud
import seaborn as sns

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


def gen_quality_h5py(bam, gff, output_filename, min_coverage=20, max_coverage=1000):
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
    gff = pd.read_csv(gff, header=None, sep='\t')

    positions = []
    pos_to_info = {}
    for line in gff.values:
        positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
        # Also keep track of the gene information
        pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()

    out_h = h5py.File(output_filename, 'w')
    for pos in tqdm(positions):
        reads = []
        try:
            for read in bam.fetch(pos):
                reads.append(read)
        except:
            print(pos)
        read_num = 0
        if len(reads) > min_coverage:
            for ri, read in enumerate(reads):
                try:
                    if read.query_sequence is not None:
                        read_name = f'{read.query_name}:{read.reference_start}-{read.reference_end}'
                        out_h.create_dataset(f'{pos}/{read_name}', data=np.array(read.query_qualities, dtype=np.int8))
                        read_num += 1
                        if read_num > max_coverage:
                            break
                except:
                    print(read.query_name)
    bam.close()
    out_h.close()

def gen_quality_data(bam, gff, output_filename=None, min_coverage=20, max_coverage=1000):
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
            for ri, read in enumerate(reads):
                if read.query_sequence is not None:
                    encodings.append([read.query_name, pos] + list(read.query_qualities))  # Want to flatten it so that we can get it easy
                    read_num += 1
                    if read_num > max_coverage:
                        break
    bam.close()
    df = pd.DataFrame(encodings)
    if output_filename is not None:
        df.to_csv(output_filename, index=False)
    return df


import math
def chunk_data(x, y, chunk_size=7): # Use 7 as this means we can relate to eligos2
    # https://stackoverflow.com/questions/13673060/split-string-into-strings-by-length
    chunks, chunk_size = len(x), chunk_size
    return [x[i:i+chunk_size] for i in range(0, chunks, chunk_size)], [y[i:i+chunk_size] for i in range(0, chunks, chunk_size)]


def gen_mapping_gene_read_dict_meta(bam_name: str, gff: str):
    """
    Counts how many reads map to a gene (using metatranscriptomics data). Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_name: list of paths to bam files
    gtf: path to reference (expect the transcriptome)
    info_cols: columns to keep from the gff info file (e.g. from that last column)
    ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey so selecting ['ID', 'Name'] would keep those two...

    Returns
    -------

    """
    gff = pd.read_csv(gff, header=None, sep='\t')
    ids = gff[0].values
    gene_info = gff[8].values
    # ToDo: generalise
    gff['gene_name'] = [f'{ids[i]}-{g.split(";")[-1].replace("ID=", "")}' for i, g in enumerate(gene_info)]
    gff['id'] = [f'{c}' for i, c in enumerate(gff[0].values)]

    # Finally map the gene name to the id
    gene_dict = dict(zip(gff['gene_name'].values, gff['id'].values))

    # Count the reads in a bam file
    all_reads = defaultdict(list)
    bam = pysam.AlignmentFile(bam_name, "rb")
    for gene in gene_dict:
        pos = gene_dict[gene]
        reads = []
        try:
            for read in bam.fetch(pos):
                reads.append(read.query_name)
                # We add the name of the read so remove any second mappings
                # We'll keep only the highest quality mapping...
                # Check if quality > 0...
                if read.mapping_quality > 0:
                    all_reads[read.query_name].append(
                        [pos, np.mean(read.query_qualities), read.mapping_quality, read.query_name])
        except:
            print(pos)
    # Keep only top mapping genes
    top_quality_gene_reads = defaultdict(list)
    total_multi = 0
    total_single = 0
    for read, reads in all_reads.items():
        if len(reads) > 1:
            total_multi += 1
        else:
            total_single += 1
        # Make only 1 map (change the index if you want something different)
        if len(reads) != 0:  # obvs don't keep if it's empty...
            top_quality_gene = max(reads, key=lambda item: item[1])
            top_quality_gene_reads[top_quality_gene[0]].append(top_quality_gene[-1])
    bam.close()
    return top_quality_gene_reads


def gen_mapping_gene_read_dict(bam_name: str, gff: str):
    """
    Counts how many reads map to a gene. Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_name: list of paths to bam files
    gtf: path to reference (expect the transcriptome)
    info_cols: columns to keep from the gff info file (e.g. from that last column)
    ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey so selecting ['ID', 'Name'] would keep those two...

    Returns
    -------

    """
    gff = pd.read_csv(gff, header=None, sep='\t')
    gene_info = gff[8].values
    # ToDo: generalise
    gff['gene_name'] = [g.split(';')[0].replace(f'ID=', '') for g in gene_info]
    starts = gff[3].values
    ends = gff[4].values
    gff['id'] = [f'{c}:{int(starts[i]) - 1}-{int(ends[i])}' for i, c in enumerate(gff[0].values)]
    # Finally map the gene name to the id
    gene_dict = dict(zip(gff['gene_name'].values, gff['id'].values))

    # Count the reads in a bam file
    qualities = []
    bam_gene_reads = {}
    all_reads = defaultdict(list)
    bam = pysam.AlignmentFile(bam_name, "rb")
    for gene in gene_dict:
        pos = gene_dict[gene]
        reads = []
        try:
            for read in bam.fetch(pos):
                reads.append(read.query_name)
                # We add the name of the read so remove any second mappings
                # We'll keep only the highest quality mapping...
                # Check if quality > 0...
                if read.mapping_quality > 0:
                    all_reads[read.query_name].append(
                        [pos, np.mean(read.query_qualities), read.mapping_quality, read.query_name])
                qualities.append(read.mapping_quality)
            bam_gene_reads[gene] = reads
        except:
            print(pos)
    # Keep only top mapping genes
    top_quality_gene_reads = defaultdict(list)
    total_multi = 0
    total_single = 0
    for read, reads in all_reads.items():
        if len(reads) > 1:
            total_multi += 1
        else:
            total_single += 1
        # Make only 1 map (change the index if you want something different)
        top_quality_gene = max(reads, key=lambda item: item[1])
        top_quality_gene_reads[top_quality_gene[0]].append(top_quality_gene[-1])
    return top_quality_gene_reads

def gen_training_h5py_position(bam: str, ref: str, positions: dict, output_filename: str,
                               min_coverage=20, max_coverage=1000):
    """
    Generate training data using a set of positions
    Parameters
    ----------
    bam: path to bam file
    ref: reference
    positions: dict of gene --> list of read ids
    output_filename
    min_coverage
    max_coverage

    Returns
    -------

    """
    bam = pysam.AlignmentFile(bam, "rb")
    fasta = pysam.FastaFile(ref)

    out_h = h5py.File(output_filename, 'w')
    for pos in tqdm(positions):
        reads = []
        try:
            for read in bam.fetch(pos):
                # Check if we want this read
                if read.query_name in positions[pos]:
                    reads.append(read)
        except:
            print(pos)
        read_num = 0
        if len(reads) > min_coverage:
            ref_str = fasta[pos]
            # Take a random sample if there are too many...
            if len(reads) > max_coverage:
                reads = random.sample(reads, max_coverage)
            for ri, read in enumerate(reads):
                try:
                    if read.query_sequence is not None:
                        seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                                   read.query_qualities)
                        read_name = read.query_name
                        seq = [ord(c) for c in seq]
                        out_h.create_dataset(f'{pos}/{read_name}/qual', data=np.array(qual, dtype=np.int8))
                        out_h.create_dataset(f'{pos}/{read_name}/seq', data=np.array(seq, dtype=np.int8))
                        out_h[f'{pos}/{read_name}'].attrs['info'] = read.reference_start
                        read_num += 1
                        if read_num > max_coverage:
                            break
                except:
                    print(read.query_name)
    bam.close()
    out_h.close()


def gen_training_h5py(bam, ref, gff, output_filename, min_coverage=20, max_coverage=1000):
    bam = pysam.AlignmentFile(bam, "rb")
    fasta = pysam.FastaFile(ref)
    gff = pd.read_csv(gff, header=None, sep='\t')

    positions = []
    pos_to_info = {}
    for line in gff.values:
        positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
        # Also keep track of the gene information
        pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()

    out_h = h5py.File(output_filename, 'w')
    for pos in tqdm(positions):
        reads = []
        try:
            for read in bam.fetch(pos):
                reads.append(read)
        except:
            print(pos)
        read_num = 0
        if len(reads) > min_coverage:
            ref_str = fasta[pos]

            for ri, read in enumerate(reads):
                try:
                    if read.query_sequence is not None:
                        seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                                   read.query_qualities)
                        read_name = read.query_name
                        # Make it totally align
                        #seq = "*" * read.reference_start + seq + "-" * (len(ref_str) - (read.reference_start + len(seq)))
                        #qual = ([0] * read.reference_start) + qual + ([0] * (len(ref_str) - (read.reference_start + len(seq))))
                        # Save sequence and quality for the read
                        seq = [ord(c) for c in seq]
                        out_h.create_dataset(f'{pos}/{read_name}/qual', data=np.array(qual, dtype=np.int8))
                        out_h.create_dataset(f'{pos}/{read_name}/seq', data=np.array(seq, dtype=np.int8))
                        out_h[f'{pos}/{read_name}'].attrs['info'] = read.reference_start
                        read_num += 1
                        if read_num > max_coverage:
                            break
                except:
                    print(read.query_name)
    bam.close()
    out_h.close()

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


def gen_coverage_over_genes_h5py(bam, ref, gff, output_filename, min_coverage=20, max_coverage=1000):
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
    out_h = h5py.File(output_filename, 'w')
    for pos in tqdm(positions):
        reads = []
        try:
            for read in bam.fetch(pos):
                reads.append(read)
        except:
            print(pos)
        read_num = 0
        if len(reads) > min_coverage:
            ref_str = fasta[pos]
            rows = []
            for ri, read in enumerate(reads):
                if read.query_sequence is not None:
                    try:
                        seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
                                                             read.query_qualities)
                        row = [read.query_name]
                        # Make it totally align
                        seq = "*" * read.reference_start + seq + "-" * (
                                len(ref_str) - (read.reference_start + len(seq)))
                        row += list(np.array(list(seq)))

                        rows.append(row)
                        read_num += 1
                        if read_num > max_coverage:
                            break
                    except:
                        print(read.query_name)
            rows = np.array(rows)
            try:
                sumarised = []
                for i in range(0, len(ref_str)):
                    sumarised.append(len(rows[rows[:, i] == 'A']))
                    sumarised.append(len(rows[rows[:, i] == 'T']))
                    sumarised.append(len(rows[rows[:, i] == 'G']))
                    sumarised.append(len(rows[rows[:, i] == 'C']))
                    sumarised.append(len(rows[rows[:, i] == '-']))  # Missing content (i.e. a deletion)
                out_h.create_dataset(f'{pos}', data=np.array(sumarised, dtype=np.int8))
                out_h[f'{pos}/{read.query_name}'].attrs['info'] = pos_to_info[pos]
                out_h[f'{pos}/{read.query_name}'].attrs['reads'] = len(reads)
            except:
                x = 1
    bam.close()
    out_h.close()

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
    output_filename = output_filename or f'reads_msa_{gene_location}.fasta'
    with open(output_filename, 'w') as fout:
        ref_str = ref[gene_location]
        fout.write(f'>Reference\n')
        fout.write(ref_str + '\n')
        for read in reads:
            seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str, read.query_qualities)
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
            # Fill them all with NA will later be filled with the mean value
            encoded[0, i] = None  # This will be dropped later ToDo: think of a better way to incorperate deletions...
            encoded[1, i] = None  # This will be dropped later ToDo: think of a better way to incorperate deletions...
            encoded[2, i] = None  # This will be dropped later ToDo: think of a better way to incorperate deletions...
            encoded[3, i] = None  # This will be dropped later ToDo: think of a better way to incorperate deletions...
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


def count_reads(bam_files: list, gff: str, info_cols=['Name', 'gene', 'gene_biotype'],
                ignore_multimapped_reads='quality'):
    """
    Counts how many reads map to a gene. Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_files: list of paths to bam files
    gtf: path to reference (expect the transcriptome)
    info_cols: columns to keep from the gff info file (e.g. from that last column)
    ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey so selecting ['ID', 'Name'] would keep those two...

    Returns
    -------

    """
    gff = pd.read_csv(gff, header=None, sep='\t')
    gene_info = gff[8].values
    # ToDo: generalise
    gff['gene_name'] = [g.split(';')[0].replace(f'ID=', '') for g in gene_info]
    starts = gff[3].values
    ends = gff[4].values
    gff['id'] = [f'{c}:{int(starts[i]) - 1}-{int(ends[i])}' for i, c in enumerate(gff[0].values)]
    # Finally map the gene name to the id
    gene_dict = dict(zip(gff['gene_name'].values, gff['id'].values))
    bam_counts = {}

    # Count the reads in a bam file
    qualities = []
    for bam_name in bam_files:
        bam_gene_reads = {}
        all_reads = defaultdict(list)
        bam = pysam.AlignmentFile(bam_name, "rb")
        for gene in gene_dict:
            pos = gene_dict[gene]
            reads = []
            try:
                for read in bam.fetch(pos):
                    reads.append(read.query_name)
                    if ignore_multimapped_reads:
                        # We add the name of the read so remove any second mappings
                        # We'll keep only the highest quality mapping...
                        # Check if quality > 0...
                        if read.mapping_quality > 0:
                            all_reads[read.query_name].append([gene, np.mean(read.query_qualities), read.mapping_quality])
                        qualities.append(read.mapping_quality)
                bam_gene_reads[gene] = reads
            except:
                print(pos)
        bam_gene_counts = {}
        # Now check for mulitmapping

        if ignore_multimapped_reads:
            top_quality_gene_reads = defaultdict(list)
            total_multi = 0
            total_single = 0
            for read, reads in all_reads.items():
                if len(reads) > 1:
                    total_multi += 1
                else:
                    total_single += 1
                # Make only 1 map (change the index if you want something different)
                top_quality_gene = max(reads, key=lambda item: item[1])
                top_quality_gene_reads[top_quality_gene[0]].append(top_quality_gene[1])
            # Now get the counts for genes
            for gene in gene_dict:
                bam_gene_counts[gene] = len(top_quality_gene_reads[gene])
            # for gene in gene_dict:
            #     # Get only reads that map a single time
            #     multimapped_reads = [r for r in bam_gene_reads[gene] if len(all_reads[r]) > 1]
            #     single_reads = [r for r in bam_gene_reads[gene] if len(all_reads[r]) == 1]
            #     # Top quality reads
            #     top_reads = [max(r, key=lambda item: item[1]) for r in bam_gene_reads[gene])
            #     bam_gene_counts[gene] = len(single_reads)
            #     total_single += len(single_reads)
            #     total_multi += len(multimapped_reads)
            print(bam_name, total_multi, total_single, total_single/(total_single + total_multi))
        else:
            for gene in gene_dict:
                bam_gene_counts[gene] = len(bam_gene_reads[gene])

        bam_counts[bam_name] = bam_gene_counts
    # Now make a df
    print(qualities)
    rows = []
    gene_to_gene_info = {}
    for info in gene_info:
        info = info.split(';')
        vals = []
        for g in info_cols:
            found = False
            for i in info:
                if g == i.split('=')[0]:
                    vals.append(i.split('=')[-1])
                    found = True
            if not found:
                vals.append(None)
        gene_to_gene_info[info[0].replace(f'ID=', '')] = vals
    for gene in gene_dict:
        row = [gene] + gene_to_gene_info[gene]
        for bam_name in bam_files:
            row.append(bam_counts[bam_name].get(gene))
        rows.append(row)
    df = pd.DataFrame(rows, columns=['gene'] + info_cols + [b.split('/')[-1].split('.')[0] for b in bam_files])
    return df

# Vis is sometimes easier to do not on the server so I'm going to look at an ORA
# of the GO terms
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
# from sciutil import SciUtil, SciException
# from sciviso import Emapplot
# import pandas as pd


def plot_clusters_go_enrichment(filename, gene_ratio_min=1, padj_max=0.05, title='GO', fig_dir='',
                    label_font_size=10, figsize=(5, 5), axis_font_size=8,save_fig=True,
                    padj_col='p.adj'):
    odds_ratio_df = pd.read_csv(filename)
    r_df = odds_ratio_df[odds_ratio_df['genes with GO and in cluster'] > gene_ratio_min]
    r_df = r_df[r_df[padj_col] < padj_max]
    r = title
    if len(r_df) > 1:
        eplot = Emapplot(r_df,
                         size_column='genes with GO and in cluster',
                         color_column='genes in cluster but not GO',
                         id_column='GO',
                         label_column='Label',
                         overlap_column='gene_names', overlap_sep=' ', title=r,
                         config={'figsize': figsize, 'label_font_size': label_font_size,
                                 'axis_font_size': axis_font_size})
        eplot.build_graph()
        plt.title(title)
        plt.gca().set_clip_on = False
        if save_fig:
            plt.savefig(f'{fig_dir}GO_{title.replace(" ", "-")}_network.svg', bbox_inches='tight',
                        transparent=True)
        plt.show()

        x, y = np.ogrid[:300, :300]

        mask = (x - 150) ** 2 + (y - 150) ** 2 > 130 ** 2
        mask = 255 * mask.astype(int)
        wordfeqs = defaultdict(int)
        for g in r_df['gene_names'].values:
            for w in g.split(' '):
                w = w.replace(' ', '.')
                wordfeqs[w] += 1
        total_words = len(wordfeqs)
        for w in wordfeqs:
            wordfeqs[w] = wordfeqs[w] / total_words
        wordcloud = WordCloud(background_color="white", mask=mask, colormap='viridis',
                              repeat=False).generate_from_frequencies(wordfeqs)

        plt.figure()
        plt.rcParams['svg.fonttype'] = 'none'  # Ensure text is saved as text
        plt.rcParams['figure.figsize'] = figsize
        font_family = 'sans-serif'
        font = 'Arial'
        sns.set(rc={'figure.figsize': figsize, 'font.family': font_family,
                    'font.sans-serif': font, 'font.size': 12}, style='ticks')
        plt.figure()
        plt.imshow(wordcloud, interpolation="bilinear")
        plt.axis("off")
        if save_fig:
            wordcloud_svg = wordcloud.to_svg(embed_font=True)
            f = open(f'{fig_dir}GO_{r}_WordCloud.svg', "w+")
            f.write(wordcloud_svg)
            f.close()
            plt.savefig(f'{fig_dir}GO_{r}_WordCloud.png', bbox_inches='tight')
        plt.show()


def run_go_enrichment(df_filename: str, go_term_filename: str, go_info_filename: str, df_gene_id='Name',
                      go_gene_id='gene', go_term='term', padj_col='padj_rna', logfc_col='logFC_rna',
                      logfc_cutoff=1.0, padj_cutoff=0.05, comparison='DE_analysis',
                      go_info_go_id='GO', go_info_description_id='Label'):
    """

    Parameters
    ----------
    df_filename: Results from DE anlaysis like DEseq2
    go_term_filename: GO file with two columns, gene & term
    go_info_filename: GO file with info for the go terms, should at least have GO id and GO name/label
    df_gene_id: column label in the df that has the gene ID in it
    go_gene_id: column label in the go DF that has the gene ID in it (needs to match)
    go_term: column label in the go DF that has the GO term in it
    padj_col: column in the DE df that has the padj in it
    logfc_col: column in the DE df that has the logFC in it
    logfc_cutoff: cutoff for logFC in the DE DF
    padj_cutoff: cutoff for padj in the DE DF
    comparison: label of the  run.
    go_info_go_id: column in the go_info_filename that has the GO id value i.e. GO:121212
    go_info_description_id: column in go_info_filename that has the description/name/label of the GO term

    Returns
    -------

    """
    df = pd.read_csv(df_filename)
    gene_go_df = pd.read_csv(go_term_filename)
    go_info_df = pd.read_csv(go_info_filename)
    columns = ['GO', 'p-value', 'odds-ratio',
               'genes with GO and in cluster',
               'genes in cluster but not GO',
               'genes in GO but not cluster',
               'genes not in GO or cluster',
               'gene_names']
    rows = []
    all_bg_genes = list(set(df[df_gene_id].values))
    term_to_label = dict(zip(go_info_df[go_info_go_id].values, go_info_df[go_info_description_id].values))
    sig = df[df[padj_col] < padj_cutoff]
    print(len(sig), len(df))
    if logfc_cutoff > 0:
        sig = sig[sig[logfc_col] > logfc_cutoff]
    else:
        sig = sig[sig[logfc_col] < logfc_cutoff]
    print(len(sig), len(df))
    cluster_genes = list(set(sig[df_gene_id].values))
    gene_id = go_gene_id
    for go in tqdm(gene_go_df[go_term].unique()):
        go_df = gene_go_df[gene_go_df[go_term] == go]
        go_enriched_genes = list(set(list(go_df[gene_id].values)))
        go_enriched_genes = [g for g in go_enriched_genes if g in all_bg_genes]
        cont_table = np.zeros((2, 2))  # Rows=ct, not ct; Cols=Module, not Modul
        if len(go_enriched_genes) > 3:  # Some DE genes
            in_go_and_cluster = len(set(cluster_genes) & set(go_enriched_genes))
            cluster_not_go = len(cluster_genes) - in_go_and_cluster
            go_not_cluster = len(go_enriched_genes) - in_go_and_cluster
            bg_genes = [g for g in all_bg_genes if g not in cluster_genes]
            if in_go_and_cluster > 3:  # Require there to be at least 3 genes overlapping
                not_go_not_cluster = len(bg_genes) - (in_go_and_cluster + cluster_not_go + go_not_cluster)
                # Populating cont table
                cont_table[0, 0] = in_go_and_cluster
                cont_table[1, 0] = cluster_not_go
                cont_table[0, 1] = go_not_cluster
                cont_table[1, 1] = not_go_not_cluster
                # Doing FET, Enrichment IN GO only.
                odds_ratio, pval = fisher_exact(cont_table, alternative="greater")
                genes_ids = list(set(cluster_genes) & set(go_enriched_genes))
                rows.append([go, pval, odds_ratio, in_go_and_cluster,
                             cluster_not_go, go_not_cluster, not_go_not_cluster,
                             ' '.join(genes_ids)])
    odds_ratio_df = pd.DataFrame(data=rows, columns=columns)
    reg, padj, a, b = multipletests(odds_ratio_df['p-value'].values,
                                    alpha=0.05, method='fdr_bh', returnsorted=False)
    odds_ratio_df['p.adj'] = padj
    ## For each gene pull out the sequence from the transcriptome file for each of our references
    odds_ratio_df['Label'] = [term_to_label.get(g) for g in odds_ratio_df['GO'].values]
    odds_ratio_df.to_csv(f'{comparison}_odds_ratio.csv', index=False)
    return odds_ratio_df


def convert_gaf_to_csv(gaf_file: str, go_term_output_file: str, go_info_output_file: str):
    """
    Converts a gaf file format to csv with term to gene so that we have easy access to these relationships

    Parameters
    ----------
    gaf_file
    go_term_output_file
    go_info_output_file

    Returns
    -------

    """

    gene_go = []
    gene_rows = []
    with open(gaf_file, 'r') as f:
        for line in f:
            if line[0] != '#' and line[0] != '!':
                line = line.split('\t')
                gene_go.append([line[2], line[4]])
                gene_rows.append(line)
    gene_go_df = pd.DataFrame(gene_go, columns=['gene', 'term'])

    all_go_df = pd.DataFrame(gene_rows,
                             columns=['UniprotID', 'ID', 'gene', 'Term thing', 'GO', 'PMID', 'Acc', 'Long', 'L',
                                      'Label', 'genes', 'Type', 'Taxon', 'taxID', 'Source', 'NA', 'NA'])

    # Remove dups incase of different sources....
    all_go_df = all_go_df.drop_duplicates()
    all_go_df.to_csv(go_info_output_file, index=False)

    gene_go_df = gene_go_df.drop_duplicates()
    gene_go_df.to_csv(go_term_output_file, index=False)


def count_reads_meta(bam_files: list, gff: str, ignore_multimapped_reads='quality', transcriptome=True):
    """
    Counts how many reads map to a gene. Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_files: list of paths to bam files
    gtf: path to reference (expect the transcriptome)
    info_cols: columns to keep from the gff info file (e.g. from that last column)
    ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey so selecting ['ID', 'Name'] would keep those two...

    Returns
    -------

    """
    gff = pd.read_csv(gff, header=None, sep='\t')
    ids = gff[0].values
    gene_info = gff[8].values
    # ToDo: generalise
    gff['gene_name'] = [f'{ids[i]}-{g.split(";")[-1].replace("ID=", "")}' for i, g in enumerate(gene_info)]
    starts = gff[3].values
    ends = gff[4].values
    if transcriptome:
        gff['id'] = [f'{c}:{int(starts[i]) - 1}-{int(ends[i])}' for i, c in enumerate(gff[0].values)]
    else:
        gff['id'] = [f'{c}' for i, c in enumerate(gff[0].values)]

    # Finally map the gene name to the id
    gene_dict = dict(zip(gff['gene_name'].values, gff['id'].values))
    bam_counts = {}

    # Count the reads in a bam file
    qualities = []
    for bam_name in bam_files:
        bam_gene_reads = {}
        all_reads = defaultdict(list)
        bam = pysam.AlignmentFile(bam_name, "rb")
        for gene in gene_dict:
            pos = gene_dict[gene]
            reads = []
            try:
                for read in bam.fetch(pos):
                    reads.append(read.query_name)
                    if ignore_multimapped_reads:
                        # We add the name of the read so remove any second mappings
                        # We'll keep only the highest quality mapping...
                        # Check if quality > 0...
                        if read.mapping_quality > 0:
                            all_reads[read.query_name].append([gene, np.mean(read.query_qualities), read.mapping_quality])
                        qualities.append(read.mapping_quality)
                bam_gene_reads[gene] = reads
            except:
                print(pos)
        bam_gene_counts = {}
        # Now check for mulitmapping
        if ignore_multimapped_reads:
            top_quality_gene_reads = defaultdict(list)
            total_multi = 0
            total_single = 0
            for read, reads in all_reads.items():
                if len(reads) > 1:
                    total_multi += 1
                else:
                    total_single += 1
                # Make only 1 map
                top_quality_gene = max(reads, key=lambda item: item[1])
                top_quality_gene_reads[top_quality_gene[0]].append(top_quality_gene[1])
            # Now get the counts for genes
            for gene in gene_dict:
                bam_gene_counts[gene] = len(top_quality_gene_reads[gene])
            try:
                print(bam_name, total_multi, total_single, total_single/(total_single + total_multi))
            except:
                print(bam_name, total_multi, total_single)
        else:
            for gene in gene_dict:
                bam_gene_counts[gene] = len(bam_gene_reads[gene])

        bam_counts[bam_name] = bam_gene_counts
    # Now make a df
    rows = []
    for gene in gene_dict:
        row = [gene]
        for bam_name in bam_files:
            row.append(bam_counts[bam_name].get(gene))
        rows.append(row)
    df = pd.DataFrame(rows, columns=['gene'] + [b.split('/')[-1].split('.')[0] for b in bam_files])
    return df


def gen_kmer_sliding_window_ref(ref: str, kmer_len: int, output_filename: str, genes=None):
    """
    Generates a positional sliding window h5 file that contains all the instances of a particular kmer.
    Get the kmers from the reference and map them to a numpy array for fast access
    -----------------------------------------------
       KMER to HDF5 here we want to iterate through the
       reference file and find where each kmer occurs
       this is then saved as kmer-->transcript-->[list of positions of kmer start]

    Parameters
    ----------
    ref: path to the reference transcriptome.
    kmer_len: length of kmer
    output_filename: string to the path of the output h5 file

    Returns
    -------

    """
    ts_dict = dict()
    # Keeps track of which transcript this kmer is in and the ID of the position it occured at (from the first position)
    with open(ref, 'r') as f:
        i = 0
        for line in f:
            if line.startswith('>'):
                ts = line.split(' ')[0][1:].strip()
                if genes is not None:
                    if ts in genes:
                        ts_dict[ts] = None
                        i += 1
                else:
                    ts_dict[ts] = None
                    i += 1
            else:
                if genes is not None:
                    if ts in genes:
                        line = line.strip()
                        ts_dict[ts] = line
                else:
                    line = line.strip()
                    ts_dict[ts] = line

    kmer_mapping = defaultdict(lambda: defaultdict(list))
    for ts, line in ts_dict.items():
        line = line
        for i in range(0, len(line) - kmer_len):  # Iterate through the sequence
            kmer = line[i:i+kmer_len]
            kmer_mapping[kmer][ts].append(i)

    kmer_h = h5py.File(output_filename, 'w')
    for k in kmer_mapping:
        for ts in kmer_mapping[k]:
            kmer_h.create_dataset(f'{k}/{ts}', data=np.array(kmer_mapping[k][ts]))

    kmer_h.close()
