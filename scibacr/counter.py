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

"""
Simple code for counting reads that mapped to genes in mapping file.
"""
import pandas as pd
import numpy as np
import pysam
from collections import defaultdict


def gen_mapping_gene_read_dict(bam_name: str, gff: str):
    """
    Counts how many reads map to a gene. Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_name: list of paths to bam files
    gff: path to reference (expect the transcriptome)

    Returns
    -------

    """
    gff = pd.read_csv(gff, header=None, sep='\t', comment="#")
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


def count_reads(bam_files: list, gff: str, info_cols=None, ignore_multimapped_reads='quality'):
    """
    Counts how many reads map to a gene. Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_files: list of paths to bam files
    gff: path to reference (expect the transcriptome)
    info_cols: columns to keep from the gff info file (e.g. from that last column) ['Name', 'gene', 'gene_biotype']
    ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey so selecting ['ID', 'Name'] would keep those two...
    ignore_multimapped_reads: whetehr or not to ignrore mulitmapped reads

    Returns
    -------

    """
    info_cols = [] if info_cols is None else info_cols
    gff = pd.read_csv(gff, header=None, sep='\t', comment="#")
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
        row = [gene, gene_dict[gene]] + gene_to_gene_info[gene]
        for bam_name in bam_files:
            row.append(bam_counts[bam_name].get(gene))
        rows.append(row)
    df = pd.DataFrame(rows, columns=['gene', 'id'] + info_cols + [b.split('/')[-1].split('.')[0] for b in bam_files])
    return df


"""

Equivalent code for 16S consolidate the two...

"""


def count_reads_meta(bam_files: list, gff: str, ignore_multimapped_reads='quality', transcriptome=True):
    """
    Counts how many reads map to a gene. Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_files: list of paths to bam files
    gff: path to reference (expect the transcriptome)
    ignore_multimapped_reads:
    transcriptome:

    Returns
    -------

    """
    gff = pd.read_csv(gff, header=None, sep='\t')
    ids = gff[0].values
    gene_info = gff[8].values
    # ToDo: generalise
    gff['gene_name'] = [g.strip() for i, g in enumerate(gene_info)]
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
        row = [gene, gene_dict[gene]]
        for bam_name in bam_files:
            row.append(bam_counts[bam_name].get(gene))
        rows.append(row)
    df = pd.DataFrame(rows, columns=['gene_name', 'id'] + [b.split('/')[-1].split('.')[0] for b in bam_files])
    return df


def count_reads_transcriptome(bam_files: list, gff: str, info_cols=None, ignore_multimapped_reads='quality'):
    """
    Counts how many reads map to a gene. Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_files: list of paths to bam files
    gff: path to reference (expect the transcriptome)
    info_cols: columns to keep from the gff info file (e.g. from that last column) ['Name', 'gene', 'gene_biotype']
    ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey so selecting ['ID', 'Name'] would keep those two...
    ignore_multimapped_reads: whetehr or not to ignrore mulitmapped reads

    Returns
    -------

    """
    info_cols = [] if info_cols is None else info_cols
    gff = pd.read_csv(gff, header=None, sep='\t', comment="#")
    gene_info = gff[8].values
    # ToDo: generalise
    gff['gene_name'] = [g.split(';')[0].replace(f'ID=', '') for g in gene_info]
    starts = gff[3].values
    ends = gff[4].values
    gff['id'] = gff['gene_name'].values #[f'{c}:{int(starts[i])}-{int(ends[i])}' for i, c in enumerate(gff[0].values)]
    # Finally map the gene name to the id
    gene_dict = dict(zip(gff['gene_name'].values, gff['id'].values))
    bam_counts = {}
    xxxx = 0
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
                xxxx += 1
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
            #print(bam_name, total_multi, total_single, total_single/(total_single + total_multi))
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
        row = [gene, gene_dict[gene]] + gene_to_gene_info[gene]
        for bam_name in bam_files:
            row.append(bam_counts[bam_name].get(gene))
        rows.append(row)
    df = pd.DataFrame(rows, columns=['gene', 'id'] + info_cols + [b.split('/')[-1].split('.')[0] for b in bam_files])
    return df



def gen_mapping_gene_read_dict_meta(bam_name: str, gff: str):
    """
    Counts how many reads map to a gene (using metatranscriptomics data). Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_name: list of paths to bam files
    gff: path to reference (expect the transcriptome)

    Returns
    -------

    """
    gff = pd.read_csv(gff, header=None, sep='\t', comment="#")
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




def gen_mapping_gene_read_dict_transcriptome(bam_name: str, gff: str):
    """
    Counts how many reads map to a gene. Output featurecounts style.

    NZ_CP092052.1   RefSeq  gene    1006731 1007435 .   +   .   ID=gene-AH5667_RS04895;Name=AH5667_RS04895;gbkey=Gene;gene_biotype=protein_coding;locus_tag=AH5667_RS04895;old_locus_tag=AH5667_000935
    Parameters
    ----------
    bam_name: list of paths to bam files
    gff: path to reference (expect the transcriptome)

    Returns
    -------

    """
    gff = pd.read_csv(gff, header=None, sep='\t', comment="#")
    gene_info = gff[8].values
    # ToDo: generalise
    gff['gene_name'] = [g.split(';')[0].replace(f'ID=', '') for g in gene_info]
    starts = gff[3].values
    ends = gff[4].values
    gff['id'] = gff['gene_name'].values #[f'{c}:{int(starts[i]) - 1}-{int(ends[i])}' for i, c in enumerate(gff[0].values)]
    # Finally map the gene name to the id
    gene_dict = dict(zip(gff['gene_name'].values, gff['id'].values))

    # Count the reads in a bam file
    qualities = []
    bam_gene_reads = {}
    all_reads = defaultdict(list)
    bam = pysam.AlignmentFile(bam_name, "rb")
    xxx = 0
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
            xxx += 1
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

