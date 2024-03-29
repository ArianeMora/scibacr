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
Contains misc functions for editing file types...
"""
import pandas as pd


def convert_transcriptome_to_gff(fasta_file: str, output_gff_file: str):
    """

    Parameters
    ----------
    fasta_file
    output_gff_file

    Returns
    -------

    """
    with open(fasta_file, 'r') as fin:
        fasta_lines = ''
        with open(output_gff_file, 'w') as fout:
            outline = ''
            line_label = ''

            for line in fin:
                if line[0] == '>':
                    if len(fasta_lines) > 0:
                        fout.write(f'{outline}\t{1}\t{len(fasta_lines)-2}\t.\t.\t.\t{line_label}\n')
                    line = line.strip()
                    # EvRCC1521_s1    EVM gene    13206   16404   .   +   .   ID=evm.EvRCC1521_s1_g1;Name=Gene_model_EvRCC1521_s1_g1
                    line = line.split(' ')
                    outline = f'{line[0][1:]}\tFA\tgene'
                    line_label = ';'.join(line).strip()[1:]
                    fasta_lines = ''
                    #>TRINITY_DN8_c1_g1_i1 len=2323 path=[0:0-2322]
                else:
                    fasta_lines += line.strip()


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

