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
import random
import pysam
import h5py
from tqdm import tqdm

from collections import defaultdict

from scibacr.misc import alignment_from_cigar


def chunk_data(x, y, chunk_size=7):  # Use 7 as this means we can relate to eligos2
    # https://stackoverflow.com/questions/13673060/split-string-into-strings-by-length
    chunks, chunk_size = len(x), chunk_size
    return [x[i:i+chunk_size] for i in range(0, chunks, chunk_size)], [y[i:i+chunk_size] for i in range(0, chunks, chunk_size)]


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


def get_kmer_encodings(kmer: str, kmer_h5: str, training_bams: list, num_samples=1000):
    """
    Gener
    Parameters
    ----------
    kmer: the kmer of interest i..e 'AAAT'
    kmer_h5: the prebuild h5 file that has all the kmer to ts --> position mappings
    training_bams: a list of paths to bams of interest
    num_samples: maximum number of reads to sample

    Returns
    -------

    """
    kmer_len = len(kmer)
    encodings = []
    nn = 0
    kmer_h5 = h5py.File(kmer_h5, 'r')
    kmer_data = kmer_h5[kmer]

    for bam in training_bams:
        run = bam.split("/")[-1].split(".")[0]
        run_h5 = h5py.File(bam, 'r')
        for gene in tqdm(kmer_data.keys()):
            # check it exists
            try:
                reads = run_h5[gene]
                if len(reads) > num_samples:
                    reads = random.sample(reads, num_samples)
                for read in reads:
                    # Get the length of the dataset
                    seq = [chr(s) for s in run_h5[gene][read]['seq']]
                    qual = [r for r in run_h5[gene][read]['qual']]
                    start = run_h5[gene][read].attrs['info']
                    # Chunk and go through each ....
                    seq = ["*"] * start + seq
                    qual = [None] * start + qual

                    for position in kmer_data[gene]:
                        kmer_seq = seq[position:position + kmer_len]
                        kmer_qual = qual[position:position + kmer_len]
                        # Potentially remove all deleted sequences...
                        if kmer_seq != ['-'] * kmer_len and kmer_seq != ['*'] * kmer_len:
                            enc = one_hot_gencode(kmer_seq, kmer_qual)
                            encodings.append([run, ''.join(kmer_seq), gene, position, position + kmer_len] + list(
                                enc.flatten()))  # Want to flatten it so that we can get it easy
            except:
                nn += 0
        print(bam, len(encodings))
    df = pd.DataFrame(encodings)
    df = df.rename(columns={0: 'Run', 1: 'kmer', 2: 'Gene', 3: 'ID', 4: 'Start'})
    return df


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
            encoded[1, i] = None
            encoded[2, i] = None
            encoded[3, i] = None
    return encoded


def create_train_chunked_set(training_files: list, kmer_len: int, genes=None, max_read_count=10):
    """

    Parameters
    ----------
    training_files
    kmer_len
    genes
    max_read_count

    Returns
    -------

    """
    encodings = []
    nn = 0
    for training_ds in training_files:
        run_h5 = h5py.File(training_ds, 'r')
        for gene in run_h5:
            if genes is None or gene in genes:
                reads = [r for r in run_h5[gene]]
                # This also stops us from getting super over represented seqs
                if len(reads) > max_read_count:
                    reads = random.sample(reads, max_read_count)
                for read in reads:
                    try:
                        # Get the length of the dataset
                        seq = [chr(s) for s in run_h5[gene][read]['seq']]
                        qual = [r for r in run_h5[gene][read]['qual']]
                        start = run_h5[gene][read].attrs['info']
                        # Chunk and go through each ....
                        seq = ["*"] * start + seq
                        qual = [None] * start + qual
                        # Chunk the data...
                        chunked_seq, chunked_qual = chunk_data(seq, qual, kmer_len)
                        # get encoder the kmer into a one hot encoded version
                        # Now we want to one hot encode the seq but add in the qual info...
                        for i, s in enumerate(chunked_seq):
                            # Remove the start and the ends
                            if None not in chunked_qual[i]:
                                enc = one_hot_gencode(s, chunked_qual[i])
                                encodings.append([training_ds, gene, i, i*kmer_len] + list(enc.flatten()))
                    except:
                        nn += 0
        run_h5.close()
    df = pd.DataFrame(encodings)
    df = df.rename(columns={0: 'Run', 1: 'Gene', 2: 'ID', 3: 'Start'})
    return df


def gen_training_h5py(bam, ref, positions: list, output_filename, min_coverage=20, max_coverage=1000):
    """
    Generates training H5 py files for a range of positions,
    Parameters
    ----------
    bam
    ref
    positions
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
