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

from scibacr.runner import Runner


class Mapper(Runner):
    """
    This class is used to run minimap2, mainly we're interested in mapping the basecalled reads
    to the reference genome or transcriptome.

    Some tools require the bed file which we can create from the output of the paf file.

    References:
        https://github.com/lh3/minimap2
    """
    def __init__(self, input_path: str, processes: str = '40', verbose=True, debug=False, logging_file=None):
        super().__init__(input_path, processes, verbose, debug, logging_file)

    def gff_to_bed(self, filter_value='') -> None:
        """
        Since all tools have a stupid amount of deps just write a simple parser for a gff file
        NC_003210.1 RefSeq  region  1   2944528 .   +   .   ID=NC_003210.1:1..2944528;Dbxref=taxon:169963;
        """
        filemap = list(self.get_mapping().keys())
        for sample in filemap:
            bedfile = self.get_bed_file(sample)
            gff_file = self.get_gtf_file(sample)
            # BED output
            # transcript:YNL036W_mRNA 0   499 cd2b3759-e172-437d-a480-f0a0eafecf75    0   +   666 998
            # chromosome, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount,
            # blockSizes, blockStarts
            with open(bedfile, 'w+') as bed:
                with open(gff_file, 'r') as gff:
                    for line in gff:
                        if line[0] == '#':
                            continue
                        else:
                            line = line.split('\t')
                            if filter_value is not None:
                                if filter_value == line[2]:
                                    # i.e. make sure that the line contains something of interest e.g. gene
                                    bed_line = [line[0], line[3], line[4], line[-1].strip(), line[5], line[6]]
                                    bed.write('\t'.join(bed_line) + '\n')
                            else:
                                bed_line = [line[0], line[3], line[4], line[-1].strip(), line[5], line[6]]
                                bed.write('\t'.join(bed_line) + '\n')

    def fasta_to_bam(self, genome=False, remove_files=True, mapping_tool='bwa') -> None:
        """
        Convert the fasta sequence to a sam file, then to a bam, finally sort this.
        We need this for later steps in the pipeline.
        """
        filemap = list(self.get_mapping().keys())
        for sample in filemap:
            if genome:
                ref = self.get_genome_reference_path(sample)  # genome
            else:
                ref = self.get_reference_path(sample)  # transcriptome
            fastq = self.get_fastq(sample)
            fasta = self.get_reference_path(sample) + '.fai'
            # Get the index (note this is expected to have already have been done)
            if genome:
                sam = f'{self.get_output_path(sample)}{sample}_genome.sam'
                bam = f'{self.get_output_path(sample)}{sample}_genome.bam'
                sorted_bam = f'{self.get_output_path(sample)}{sample}_genome.sorted.bam'
            else:
                sam = f'{self.get_output_path(sample)}{sample}.sam'
                bam = f'{self.get_output_path(sample)}{sample}.bam'
                sorted_bam = f'{self.get_output_path(sample)}{sample}.sorted.bam'
            # Run sam
            if genome:
                # ../software/bwa/./bwa mem tRNA_ref_new.fasta -W 13 -k 6 -x ont2d ecoli_tRNA/ERR6751707.fastq.gz > ERR6751707.sam
                if mapping_tool == 'bwa': # -k 6 -x ont2d
                    self.run(['../software/bwa/./bwa', 'mem',  ref, '-W', '13', '-k ', '6', '-x', 'ont2d', fastq, '>', sam])
                else:
                    self.run(['minimap2', '-ax', 'splice', '-uf', '-k14 ', '-t', '20', ref, fastq, '>', sam])
            else:  # add pre- prefix because there are sometimes issues with the size when using transcriptome
                # https://github.com/lh3/minimap2/issues/301
                # https://github.com/lh3/minimap2/issues/37
                # https://github.com/lh3/minimap2/issues/358 #  '--split-prefix=pre', didn't work?
                if mapping_tool == 'bwa':
                    self.run(['../software/bwa/./bwa', '-W', '13', '-k ', '6', '-x', 'ont2d', ref, fastq, '>', sam])
                else:
                    self.run(['minimap2', '-ax', 'splice', '-uf', '-k14 ', '-t', '20', ref, fastq, '>', sam])
            # Run bam
            self.run(['samtools', 'view', '-S', sam, '-bh', '-t', fasta, '>', bam])
            # Sort bam
            self.run(['samtools', 'sort', bam, '>', sorted_bam])
            # Index bam
            self.run(['samtools', 'index', sorted_bam])
            # Remove other files
            if remove_files:
                self.run(['rm', bam])  # Just leave the sorted bam
                self.run(['rm', sam])

    def genome_to_transcripts(self) -> None:
        """
        First step is to make sure that we have a reference transcriptome.
        """
        filemap = list(self.get_mapping().keys())
        for sample in filemap:
            ref = self.get_genome_reference_path(sample)
            # Index bam
            self.run(['bedtools', 'getfasta', '-fi', ref, '-bed',
                      self.get_gtf_file(sample), '-fo', self.get_reference_path(sample)])

