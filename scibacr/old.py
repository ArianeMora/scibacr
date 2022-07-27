
# def gen_quality_data(bam, gff, output_filename=None, min_coverage=20, max_coverage=1000):
#     """
#
#     Parameters
#     ----------
#     bam
#     gff
#     min_coverage
#     max_coverage
#
#     Returns
#     -------
#
#     """
#     bam = pysam.AlignmentFile(bam, "rb")
#     gff = pd.read_csv(gff, header=None, sep='\t')
#
#     positions = []
#     pos_to_info = {}
#     for line in gff.values:
#         positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
#         # Also keep track of the gene information
#         pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()
#
#     encodings = []
#     for pos in tqdm(positions):
#         reads = []
#         for read in bam.fetch(pos):
#             reads.append(read)
#         read_num = 0
#         if len(reads) > min_coverage:
#             for ri, read in enumerate(reads):
#                 if read.query_sequence is not None:
#                     encodings.append([read.query_name, pos] + list(read.query_qualities))  # Want to flatten it so that we can get it easy
#                     read_num += 1
#                     if read_num > max_coverage:
#                         break
#     bam.close()
#     df = pd.DataFrame(encodings)
#     if output_filename is not None:
#         df.to_csv(output_filename, index=False)
#     return df
#
#
#
# def gen_training(bam, ref, gff, output_filename=None, kmer_len=7, min_coverage=20, max_coverage=1000):
#     bam = pysam.AlignmentFile(bam, "rb")
#     fasta = pysam.FastaFile(ref)
#     gff = pd.read_csv(gff, header=None, sep='\t')
#
#     positions = []
#     pos_to_info = {}
#     for line in gff.values:
#         positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
#         # Also keep track of the gene information
#         pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()
#
#     encodings = []
#
#     for pos in tqdm(positions):
#         reads = []
#         for read in bam.fetch(pos):
#             reads.append(read)
#         read_num = 0
#         if len(reads) > min_coverage:
#             ref_str = fasta[pos]
#
#             for ri, read in enumerate(reads):
#                 if read.query_sequence is not None:
#                     seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
#                                                                read.query_qualities)
#                     name = read.query_name
#                     # Make it totally align
#                     seq = "*" * read.reference_start + seq + "-" * (len(ref_str) - (read.reference_start + len(seq)))
#                     qual = ([None] * read.reference_start) + qual + (
#                                 [None] * (len(ref_str) - (read.reference_start + len(seq))))
#                     # Chunk the data...
#                     chunked_seq, chunked_qual = chunk_data(seq, qual, kmer_len)
#
#                     # Now we want to one hot encode the seq but add in the qual info...
#                     for i, s in enumerate(chunked_seq):
#                         enc = one_hot_gencode(s, chunked_qual[i])
#                         encodings.append(
#                             [name, pos, i * kmer_len] + list(
#                                 enc.flatten()))  # Want to flatten it so that we can get it easy
#                     read_num += 1
#                     if read_num > max_coverage:
#                         break
#     bam.close()
#     df = pd.DataFrame(encodings)
#     df = df.dropna()
#     if output_filename is not None:
#         df.to_csv(output_filename, index=False)
#     return df
#
#
# def gen_quality_h5py(bam, gff, output_filename, min_coverage=20, max_coverage=1000):
#     """
#
#     Parameters
#     ----------
#     bam
#     gff
#     min_coverage
#     max_coverage
#
#     Returns
#     -------
#
#     """
#     bam = pysam.AlignmentFile(bam, "rb")
#     gff = pd.read_csv(gff, header=None, sep='\t')
#
#     positions = []
#     pos_to_info = {}
#     for line in gff.values:
#         positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
#         # Also keep track of the gene information
#         pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()
#
#     out_h = h5py.File(output_filename, 'w')
#     for pos in tqdm(positions):
#         reads = []
#         try:
#             for read in bam.fetch(pos):
#                 reads.append(read)
#         except:
#             print(pos)
#         read_num = 0
#         if len(reads) > min_coverage:
#             for ri, read in enumerate(reads):
#                 try:
#                     if read.query_sequence is not None:
#                         read_name = f'{read.query_name}:{read.reference_start}-{read.reference_end}'
#                         out_h.create_dataset(f'{pos}/{read_name}', data=np.array(read.query_qualities, dtype=np.int8))
#                         read_num += 1
#                         if read_num > max_coverage:
#                             break
#                 except:
#                     print(read.query_name)
#     bam.close()
#     out_h.close()
#
#
#
# def gen_coverage_over_genes_h5py(bam, ref, gff, output_filename, min_coverage=20, max_coverage=1000):
#     """
#
#     Parameters
#     ----------
#     bam
#     ref
#     gff
#     min_coverage
#     max_coverage
#
#     Returns
#     -------
#
#     """
#     bam = pysam.AlignmentFile(bam, "rb")
#     fasta = pysam.FastaFile(ref)
#     gff = pd.read_csv(gff, header=None, sep='\t')
#
#     positions = []
#     pos_to_info = {}
#     for line in gff.values:
#         positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
#         # Also keep track of the gene information
#         pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()
#     out_h = h5py.File(output_filename, 'w')
#     for pos in tqdm(positions):
#         reads = []
#         try:
#             for read in bam.fetch(pos):
#                 reads.append(read)
#         except:
#             print(pos)
#         read_num = 0
#         if len(reads) > min_coverage:
#             ref_str = fasta[pos]
#             rows = []
#             for ri, read in enumerate(reads):
#                 if read.query_sequence is not None:
#                     try:
#                         seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
#                                                              read.query_qualities)
#                         row = [read.query_name]
#                         # Make it totally align
#                         seq = "*" * read.reference_start + seq + "-" * (
#                                 len(ref_str) - (read.reference_start + len(seq)))
#                         row += list(np.array(list(seq)))
#
#                         rows.append(row)
#                         read_num += 1
#                         if read_num > max_coverage:
#                             break
#                     except:
#                         print(read.query_name)
#             rows = np.array(rows)
#             try:
#                 sumarised = []
#                 for i in range(0, len(ref_str)):
#                     sumarised.append(len(rows[rows[:, i] == 'A']))
#                     sumarised.append(len(rows[rows[:, i] == 'T']))
#                     sumarised.append(len(rows[rows[:, i] == 'G']))
#                     sumarised.append(len(rows[rows[:, i] == 'C']))
#                     sumarised.append(len(rows[rows[:, i] == '-']))  # Missing content (i.e. a deletion)
#                 out_h.create_dataset(f'{pos}', data=np.array(sumarised, dtype=np.int8))
#                 out_h[f'{pos}'].attrs['info'] = pos_to_info[pos]
#                 out_h[f'{pos}'].attrs['reads'] = len(reads)
#             except:
#                 x = 1
#     bam.close()
#     out_h.close()
#
#
#
# def gen_coverage_over_genes(bam, ref, gff, min_coverage=20, max_coverage=1000):
#     """
#
#     Parameters
#     ----------
#     bam
#     ref
#     gff
#     min_coverage
#     max_coverage
#
#     Returns
#     -------
#
#     """
#     bam = pysam.AlignmentFile(bam, "rb")
#     fasta = pysam.FastaFile(ref)
#     gff = pd.read_csv(gff, header=None, sep='\t')
#
#     positions = []
#     pos_to_info = {}
#     for line in gff.values:
#         positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
#         # Also keep track of the gene information
#         pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()
#
#     gene_rows = []
#     for pos in tqdm(positions):
#         reads = []
#         for read in bam.fetch(pos):
#             reads.append(read)
#         read_num = 0
#         if len(reads) > min_coverage:
#             ref_str = fasta[pos]
#             rows = []
#             for ri, read in enumerate(reads):
#                 if read.query_sequence is not None:
#                     seq, qual, ref, ins = alignment_from_cigar(read.cigartuples, read.query_sequence,
#                                                                ref_str, read.query_qualities)
#                     row = [read.query_name]
#                     # Make it totally align
#                     seq = "*" * read.reference_start + seq + "-" * (
#                             len(ref_str) - (read.reference_start + len(seq)))
#                     row += list(np.array(list(seq)))
#                     rows.append(row)
#                     read_num += 1
#                     if read_num > max_coverage:
#                         break
#             rows = np.array(rows)
#             try:
#                 sumarised = []  # np.zeros(len(ref_str)*4)
#                 for i in range(0, len(ref_str)):
#                     sumarised.append(len(rows[rows[:, i] == 'A']))
#                     sumarised.append(len(rows[rows[:, i] == 'T']))
#                     sumarised.append(len(rows[rows[:, i] == 'G']))
#                     sumarised.append(len(rows[rows[:, i] == 'C']))
#                     sumarised.append(len(rows[rows[:, i] == '-']))  # Missing content (i.e. a deletion)
#                 sumarised_arr = np.array(sumarised)
#                 gene_rows.append([pos, pos_to_info[pos], np.mean(sumarised_arr), np.max(sumarised_arr),
#                                   np.var(sumarised_arr)] + sumarised)
#             except:
#                 x = 1
#     df = pd.DataFrame(gene_rows)  # Can't give it columns since there will be differing column lengths
#     bam.close()
#     return df
#
# def alignment_from_cigar_inc_inserts(cigar: str, alignment: str, ref: str):
#     """
#     Generate the alignment from the cigar string.
#     Operation	Description	Consumes query	Consumes reference
#     0 M	alignment match (can be a sequence match or mismatch)	yes	yes
#     1 I	insertion to the reference	yes	no
#     2 D	deletion from the reference	no	yes
#     3 N	skipped region from the reference	no	yes
#     4 S	soft clipping (clipped sequences present in SEQ)	yes	no
#     5 H	hard clipping (clipped sequences NOT present in SEQ)	no	no
#     6 P	padding (silent deletion from padded reference)	no	no
#     7 =	sequence match	yes	yes
#     8 X	sequence mismatch	yes	yes
#
#     Parameters
#     ----------
#     cigar: 49S9M1I24M2D9M2I23M1I3M1D13M1D10M1D13M3D12M1D23M1
#     alignment: GCTGATCACAACGAGAGCTCTCGTTGCTCATTACCCCTAAGGAACTCAAATGACGGTTAAAAACTTGTTTTGCT
#     ref: reference string (as above but from reference)
#
#     Returns
#     -------
#
#     """
#     new_seq = ''
#     ref_seq = ''
#     pos = 0
#     ref_pos = 0
#     for op, op_len in cigar:
#         if op == 0:  # alignment match (can be a sequence match or mismatch)
#             new_seq += alignment[pos:pos+op_len]
#             ref_seq += ref[ref_pos:ref_pos + op_len]
#             pos += op_len
#             ref_pos += op_len
#         elif op == 1:  # insertion to the reference
#             new_seq += alignment[pos:pos + op_len]
#             ref_seq += '^'*op_len
#             pos += op_len
#         elif op == 2:  # deletion from the reference
#             new_seq += '-'*op_len
#             ref_seq += ref[ref_pos:ref_pos + op_len]
#             ref_pos += op_len
#         elif op == 3:  # skipped region from the reference
#             new_seq += '*'*op_len
#             ref_pos += op_len
#         elif op == 4:  # soft clipping (clipped sequences present in SEQ)
#             new_seq += alignment[pos:pos + op_len]
#             ref_seq += '>' * op_len
#             pos += op_len
#         elif op == 5:  # hard clipping (clipped sequences NOT present in SEQ)
#             continue
#         elif op == 6:  # padding (silent deletion from padded reference)
#             continue
#         elif op == 7:  # sequence mismatch
#             new_seq += alignment[pos:pos + op_len]
#             ref_seq += ref[ref_pos:ref_pos + op_len]
#             pos += op_len
#             ref_pos += op_len
#     return new_seq, ref_seq
#
#
#
# def gen_training_h5py(bam, ref, gff, output_filename, min_coverage=20, max_coverage=1000):
#     bam = pysam.AlignmentFile(bam, "rb")
#     fasta = pysam.FastaFile(ref)
#     gff = pd.read_csv(gff, header=None, sep='\t')
#
#     positions = []
#     pos_to_info = {}
#     for line in gff.values:
#         positions.append(f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}')
#         # Also keep track of the gene information
#         pos_to_info[f'{line[0]}:{int(line[3]) - 1}-{int(line[4])}'] = line[8].strip()
#
#     out_h = h5py.File(output_filename, 'w')
#     for pos in tqdm(positions):
#         reads = []
#         try:
#             for read in bam.fetch(pos):
#                 reads.append(read)
#         except:
#             print(pos)
#         read_num = 0
#         if len(reads) > min_coverage:
#             ref_str = fasta[pos]
#
#             for ri, read in enumerate(reads):
#                 try:
#                     if read.query_sequence is not None:
#                         seq, ref, qual, ins = alignment_from_cigar(read.cigartuples, read.query_sequence, ref_str,
#                                                                    read.query_qualities)
#                         read_name = read.query_name
#                         # Save sequence and quality for the read
#                         seq = [ord(c) for c in seq]
#                         out_h.create_dataset(f'{pos}/{read_name}/qual', data=np.array(qual, dtype=np.int8))
#                         out_h.create_dataset(f'{pos}/{read_name}/seq', data=np.array(seq, dtype=np.int8))
#                         out_h[f'{pos}/{read_name}'].attrs['info'] = read.reference_start
#                         read_num += 1
#                         if read_num > max_coverage:
#                             break
#                 except:
#                     print(read.query_name)
#     bam.close()
#     out_h.close()