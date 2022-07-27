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

import subprocess
import pandas as pd
import time
import os


class Runner:
    """
    Base runner class that has the default getter and setters and also the default setup methods.

    Here we can simply change the running format & logging for all the other 'Runners'.
    """

    def __init__(self, input_path: str, processes: str = '40', verbose=True, debug=False, logging_file=None):

        self.info = None
        self._processes = processes
        self._verbose = verbose
        self._debug = debug
        self._logging_file = logging_file
        self.df = None
        self.input_path = input_path
        self.name = 'Runner'
        self.setup(input_path)

    def save(self, output_path: str):
        """
        Here we just
        """
        output_path = output_path if output_path else ""
        # Convert the dictionary back to a DF and save with all the new columns
        updated_df = pd.DataFrame.from_dict(self.info, orient='index')
        updated_df.to_csv(f'{output_path}{self.name}.csv', sep=',', index=False)

    def setup(self, input_path: str):
        """
        Setup the dictionary based on the users input for easy access.
        We also create the names of the output files here. They will be accessible by all runners.
        """

        df = pd.read_csv(input_path)
        basecalls = df['basecall_fastq_path'].tolist()
        reference = df['reference_fasta_path'].tolist()
        output_path = df['output_folder'].tolist()
        gtf_file = df['gtf_file'].tolist()

        ds_dict = {}
        for i, sample in enumerate(df['Sample_name']):
            out = output_path[i] if output_path[i] else ''
            ds_dict[sample] = {'basecall_fastq_path': basecalls[i],
                               'genome_reference_fasta_path': reference[i],
                               'transcript_reference_fasta_path': ''.join(reference[i].replace(".fna", "").replace(".fasta", "")) + '_transcripts.fasta',
                               'output_folder': out,
                               'bedfile': f'{out}{sample}.bed',
                               'paffile': f'{out}{sample}.paf',
                               'transcript_bam': f'{out}{sample}.sorted.bam',
                               'genome_bam': f'{out}{sample}_genome.sorted.bam',
                               'gtf_file': gtf_file[i],
                               'Sample_name': sample}
        self.info = ds_dict
        self.df = df  # Save this and we'll update it at the end and save it to a new file

    def run(self, command):
        # Check if there is a > or | since we'll change it to a string
        # possibly hacky.
        if isinstance(command, list):
            if '>' in command or '|' in command:
                command = ' '.join(command)
        if self._debug:
            self.printv(command)
        else:
            # Keep track of the time and resource usage
            if self._logging_file:
                start = time.time()
                with open(self._logging_file, 'a') as f:
                    if isinstance(command, str):
                        os.system(command)  # More vunerable
                        end = time.time()
                        f.write(f'{self.name}\t{end - start}\t{command}\tNA\n')
                    else:
                        # Run process and get the output
                        p = str(subprocess.check_output(command)).replace('\t', ' ')
                        end = time.time()
                        f.write(f'{self.name}\t{end - start}\t{" ".join(command)}\t{p}\n')
            else:
                if isinstance(command, str):
                    os.system(command)
                elif isinstance(command, list):
                    subprocess.run(command)
                else:
                    print("MUST BE LIST OR STRING")

    def printv(self, *args):
        if self._verbose:
            print(*args)
        return

    """
    -----------------------------------------------------------------
    Simple getters and setters.
    -----------------------------------------------------------------
    """
    def get_mapping(self):
        """Return a dictionary with sample instances as the keys and filepath instances as the values."""
        return self.info

    def get_reference_path(self, sample: str) -> str:
        """This is the transcriptome. Used in Tombo and minimap etc."""
        return self.info[sample]['transcript_reference_fasta_path']

    def get_genome_reference_path(self, sample: str) -> str:
        """This is the genome, used in Eligos."""
        return self.info[sample]['genome_reference_fasta_path']

    def get_output_path(self, sample: str) -> str:
        """Return the sample output dir filepath instances from input dictionary."""
        return self.info[sample]['output_folder']

    def get_fastq(self, sample: str) -> str:
        """ Get the fastq file for a sample. """
        return self.info[sample]['basecall_fastq_path']

    def get_bam_path(self, sample: str) -> str:
        """Given we're interested in RNA use transcript bam for the majority."""
        return self.info[sample]['transcript_bam']

    def get_genome_mapped_bam_path(self, sample: str) -> str:
        """Again, eligos requires a mapping between genome with the bed file so use this."""
        return self.info[sample]['genome_bam']

    def get_bed_file(self, sample):
        """Return the sample bedfile dir filepath instances from updated dictionary."""
        return self.info[sample]['bedfile']

    def get_paf_file(self, sample):
        """Return the sample bedfile dir filepath instances from updated dictionary."""
        return self.info[sample]['paffile']

    def get_fast5_path(self, sample):
        """Return the path for the single read fast5 files."""
        return self.info[sample]['fast5_path']

    def get_seq_summary_path(self, sample):
        """ Return the output for the sequencing summary from guppy."""
        return self.info[sample]['sequence_summary_path']

    def get_file_mapping(self, sample):
        """ Get the filemapping from nanocompore"""
        return self.info[sample]['file_mapping']

    def get_gtf_file(self, sample):
        """ Get the filemapping from nanocompore"""
        return self.info[sample]['gtf_file']

