# m_matschiner Wed May 19 18:38:31 CEST 2021

# This script reads a MAF file and generates a set of alignments
# in Phylip format for regions sampled across the chromosome.

class Window(object):

    def __init__(self, ids, seqs, scf, start_pos, end_pos, n_polymorphic):

        # Prepare an alignment without variation from randomly chosen nucleotides.
        self.ids = ids
        self.seqs = seqs
        self.scf = scf
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.n_polymorphic = n_polymorphic

    def write_align(self, path_prefix, file_format):

        # Get the maximum id length.
        max_id_length = 0
        for idx in self.ids:
            if len(idx) > max_id_length:
                max_id_length = len(self.ids)

        # Prepare the alignment string in Phylip format.
        if file_format == "phylip":
            outstring = str(len(self.seqs)) + " " + str(len(self.seqs[0])) + "\n"
            for y in range(len(self.seqs)):
                outstring += self.ids[y].ljust(max_id_length+2) + "".join(self.seqs[y]) + "\n"
            f_name = path_prefix + "_" + self.scf + "_" + str(self.start_pos) + "_" + str(self.end_pos) + ".phy"
        elif file_format == "fasta":
            outstring = ""
            for y in range(len(self.seqs)):
                outstring += ">" + self.ids[y] + "\n"
                outstring += "".join(self.seqs[y]) + "\n"
            f_name = path_prefix + "_" + self.scf + "_" + str(self.start_pos) + "_" + str(self.end_pos) + ".fasta"
        elif file_format == "nexus":
            outstring = "#nexus\n"
            outstring += "begin data;\n"
            outstring += "  dimensions ntax=" + str(len(self.ids)) + " nchar=" + str(len(self.seqs[0])) + ";\n"
            outstring += "  format datatype=dna missing=? gap=-;\n"
            outstring += "  matrix\n"
            for y in range(len(self.seqs)):
                outstring += "  " + self.ids[y].ljust(max_id_length+2) + "".join(self.seqs[y]) + "\n"
            outstring += "  ;\n"
            outstring += "end;\n"
            f_name = path_prefix + "_" + self.scf + "_" + str(self.start_pos) + "_" + str(self.end_pos) + ".nex"
            
        # Write the output file.
        with open(f_name, "w") as f:
            f.write(outstring)
            print("Wrote file " + f_name + " with a length of " + str(len(self.seqs[0])) + " bp and " + str(self.n_polymorphic) + " polymorphic sites.")



# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse
import textwrap
from random import randint
from random import sample
import copy
from pathlib import Path

# Set up the argument parser.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    ------------------------------------------------------------
      This script reads a MAF file and generates a set of alignments
      in Phylip format for regions sampled across the chromosome.

      Run e.g. with
      python3 make_alignments_from_maf.py input.maf output_prefix -n 1000 -l 1000
    '''))
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.1'
    )
parser.add_argument(
    '-n',
    nargs=1,
    type=int,
    default=[1],
    dest='max_n_align',
    help='Maximum number of alignments (default: 1)'
        )
parser.add_argument(
    '-l',
    nargs=1,
    type=int,
    default=[1],
    dest='min_align_length',
    help='Minimum length of alignments in bp (default: 1)'
        )
parser.add_argument(
    '-k',
    nargs=1,
    type=str,
    default=["None"],
    dest='max_align_length',
    help='Maximum length of alignments in bp (default: None)'
        )
parser.add_argument(
    '-c',
    nargs=1,
    type=float,
    default=[0.0],
    dest='min_completeness',
    help='Minimum completeness of alignments as a proportion (default: 0)'
        )
parser.add_argument(
    '-x',
    nargs=1,
    type=int,
    default=[1],
    dest='min_n_seqs',
    help='Minimum number of sequences in an alignment (default: 1)'
        )
parser.add_argument(
    '-m',
    nargs=1,
    type=int,
    default=[0],
    dest='min_n_polymorphic',
    help='Minimum number of polymorphic sites in an alignment (default: 0)'
        )
parser.add_argument(
    '-p',
    nargs=1,
    type=str,
    default=["./"],
    dest='path',
    help='Path to which to write output files (default: ./)'
        )
parser.add_argument(
    '-f',
    nargs=1,
    type=str,
    default=["phylip"],
    dest='format',
    help='Format for output files (default: phylip; other options: fasta, nexus)'
        )
parser.add_argument(
    'maf',
    nargs=1,
    type=str,
    help='name of the input file in MAF format'
    )
parser.add_argument(
    'prefix',
    nargs=1,
    type=str,
    help='prefix for output file names'
    )

# Parse command-line arguments.
args = parser.parse_args()
maf_name = args.maf[0]
prefix = args.prefix[0]
max_n_align = args.max_n_align[0]
min_n_seqs = args.min_n_seqs[0]
min_align_length = args.min_align_length[0]
max_align_length = args.max_align_length[0]
if max_align_length == "None":
    max_align_length = -1
else:
    max_align_length = int(max_align_length)
min_completeness = args.min_completeness[0]
min_n_polymorphic = args.min_n_polymorphic[0]
path = args.path[0] + "/"
path = path.replace("//","/")
Path(path).mkdir(parents=True, exist_ok=True)
path_prefix = path + prefix
file_format = args.format[0].lower()

# Make sure that the file format is known.
if file_format not in ["phylip", "fasta", "nexus"]:
    print("ERROR: Unknown file format " + file_format + "!")
    sys.exit(1)

# Parse the MAF input file line by line.
windows = []
a_line = None
s_lines = []
first_set_of_ids = None
with open(maf_name) as maf_file:
    for line in maf_file:
        if line[:2] == "a\n":
            a_line = line.strip()
        elif line.strip() == "" and a_line != None:
            ids = []
            scfs = []
            start_poss = []
            end_poss = []
            strands = []
            seqs = []
            for other_line in s_lines:
                ary = other_line.split()
                ids.append(ary[1].split(".")[0]) # This assumes that IDs are in the format "species.scaffold"
                scfs.append(ary[1].split(".")[1])
                start_poss.append(int(ary[2]))
                end_poss.append(int(ary[2]) + int(ary[3]))
                strands.append(ary[4])
                seqs.append(ary[6])
            # Use only alignments with sufficient sequences.
            if len(ids) >= min_n_seqs:
                # Make sure that the ids are unique to exclude paralogs.
                if len(set(ids)) == len(ids):
                    # Make sure that all sequences have the same length.
                    for seq in seqs:
                        if len(seq) != len(seqs[0]):
                            print("ERROR: Sequences differ in length!")
                            sys.exit(1)
                    # Make sure that the ids are unchanged from the first alignment.
                    if first_set_of_ids == None:
                        first_set_of_ids = ids
                    else:
                        if ids != first_set_of_ids:
                            print("ERROR: The order of IDs seems to vary between alignments!")
                            sys.exit(0)
                    # Use only alignments with sufficient length.
                    if len(seqs[0]) >= min_align_length:
                        # Calculate the completeness of the alignment.
                        n_non_missing = 0
                        n_missing = 0
                        for seq in seqs:
                            n_non_missing += seq.count("A") + seq.count("a")
                            n_non_missing += seq.count("C") + seq.count("c")
                            n_non_missing += seq.count("G") + seq.count("g")
                            n_non_missing += seq.count("T") + seq.count("t")
                            n_missing += seq.count("-") + seq.count("N") + seq.count("n") + seq.count("?")
                        completeness = n_non_missing / (n_non_missing + n_missing)
                        # Make sure that all characters are counted.
                        if n_non_missing + n_missing != len(ids) * len(seqs[0]):
                            print("ERROR: Not all characters are counted!")
                            sys.exit(1)
                        # Use only alignments with sufficient completeness.
                        if completeness >= min_completeness:
                            window_rel_start_poss = []
                            window_rel_end_poss = []
                            # If a maximum length has been specified, extract one or more windows from it.
                            if max_align_length > 1 and len(seqs[0]) >= max_align_length:
                                n_alignments_from_this_block = int(len(seqs[0])/max_align_length)
                                window_center_spacing = int(len(seqs[0]) / n_alignments_from_this_block)
                                window_rel_start_pos = int(window_center_spacing/2) - int(max_align_length/2)
                                window_rel_end_pos = window_rel_start_pos + max_align_length
                                window_rel_start_poss.append(window_rel_start_pos)
                                window_rel_end_poss.append(window_rel_end_pos)
                                for _ in range(n_alignments_from_this_block-1):
                                    window_rel_start_pos += window_center_spacing
                                    window_rel_start_poss.append(window_rel_start_pos)
                                    window_rel_end_pos += window_center_spacing
                                    window_rel_end_poss.append(window_rel_end_pos)
                            else:
                                window_rel_start_poss = [0]
                                window_rel_end_poss = [len(seqs[0])]
                            # For each window, truncate the sequences, and check if they're still polymorphic enough, then write the alignment.
                            for x in range(len(window_rel_start_poss)):
                                window_seqs = []
                                for seq in seqs:
                                    window_seqs.append(seq[window_rel_start_poss[x]:window_rel_end_poss[x]])
                                # Calculate the number of polymorphic sites in the alignment.
                                n_polymorphic = 0
                                for pos in range(len(window_seqs[0])):
                                    chars_this_site = []
                                    for window_seq in window_seqs:
                                        if window_seq[pos] in ["A", "a", "C", "c", "G", "g", "T", "t"]:
                                            chars_this_site.append(window_seq[pos].upper())
                                    if len(set(chars_this_site)) > 1:
                                        n_polymorphic += 1
                                # Use only alignments with sufficient numbers of polymorphic sites:
                                if n_polymorphic >= min_n_polymorphic:
                                    windows.append(Window(ids, window_seqs, scfs[0], start_poss[0] + window_rel_start_poss[x], end_poss[0] + window_rel_end_poss[x], n_polymorphic))
            a_line = None
            s_lines = []
        elif line[0:2] == "s\t" and a_line != None:
            s_lines.append(line.strip())

# Write alignments to phylip files, considering the specified maximum number of alignments.
n_aligns_written = 0
if len(windows) <= max_n_align:
    for window in windows:
        window.write_align(path_prefix, file_format)
        n_aligns_written += 1
else:   
    for window in sample(windows, max_n_align):
        window.write_align(path_prefix, file_format)
        n_aligns_written += 1
print("Wrote " + str(n_aligns_written) + " out of " + str(len(windows)) + " alignments to directory " + path + ".")
