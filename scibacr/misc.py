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


def alignment_from_cigar(cigar: str, alignment: str, ref: str) -> tuple[str, str, list]:
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
    inserts = []
    pos = 0
    ref_pos = 0
    for op, op_len in cigar:
        if op == 0:  # alignment match (can be a sequence match or mismatch)
            new_seq += alignment[pos:pos+op_len]
            ref_seq += ref[ref_pos:ref_pos + op_len]
            pos += op_len
            ref_pos += op_len
        elif op == 1:  # insertion to the reference
            inserts.append(alignment[pos - 1:pos+op_len])
            pos += op_len
        elif op == 2:  # deletion from the reference
            new_seq += '-'*op_len
            ref_seq += ref[ref_pos:ref_pos + op_len]
            ref_pos += op_len
        elif op == 3:  # skipped region from the reference
            new_seq += '*'*op_len
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
            pos += op_len
            ref_pos += op_len
    return new_seq, ref_seq, inserts


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