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
import random

import numpy as np
from tqdm import tqdm
from scibacr.misc import alignment_from_cigar
import matplotlib.pyplot as plt
from collections import defaultdict
from wordcloud import WordCloud
import seaborn as sns

# Vis is sometimes easier to do not on the server so I'm going to look at an ORA
# of the GO terms
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sciviso import Emapplot
import pandas as pd

"""
Functions for visualisation and interpretation of data
"""


def write_msa_over_gene(gene_location: str, bam, ref, output_filename=None, max_reads=200, read_dict=None):
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
    if read_dict is not None:
        reads = []
        for read in bam.fetch(gene_location):
            # Check if we want this read
            if read.query_name in read_dict[gene_location]:
                reads.append(read)
    else:
        reads = []
        for read in bam.fetch(gene_location):
            reads.append(read)
    if len(reads) > max_reads:
        reads = random.sample(reads, max_reads)
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
