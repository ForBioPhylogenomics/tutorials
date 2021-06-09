# m_matschiner Wed May 19 18:38:31 CEST 2021

# This script reads a MAF file and generates a file in VCF
# format.

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
from datetime import date

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
    '-l',
    nargs=1,
    type=int,
    default=[1],
    dest='min_align_length',
    help='Minimum length of alignment blocks in bp (default: 1)'
        )
parser.add_argument(
    '-c',
    nargs=1,
    type=float,
    default=[0.0],
    dest='min_completeness',
    help='Minimum completeness of alignment blocks as a proportion (default: 0)'
        )
parser.add_argument(
    '-x',
    nargs=1,
    type=int,
    default=[1],
    dest='min_n_seqs',
    help='Minimum number of sequences in alignment blocks (default: 1)'
        )
parser.add_argument(
    'maf',
    nargs=1,
    type=str,
    help='name of the input file in MAF format'
    )
parser.add_argument(
    'vcf',
    nargs=1,
    type=str,
    help='name of the output file in VCF format'
    )

# Parse command-line arguments.
args = parser.parse_args()
maf_name = args.maf[0]
vcf_name = args.vcf[0]
min_n_seqs = args.min_n_seqs[0]
min_align_length = args.min_align_length[0]
min_completeness = args.min_completeness[0]
date = date.today()

# Parse the MAF input file line by line.
a_line = None
s_lines = []
first_set_of_ids = None
vcf_string = "##fileformat=VCFv4.2\n"
vcf_string += "##" + str(date.year) + str(date.month) + str(date.day) + "\n"
vcf_string += "##source=make_vcf_from_maf.py\n"
vcf_string += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
n_variants = 0
with open(maf_name) as maf_file:
    for line in maf_file:
        if line[:2] == "a\n":
            a_line = line.strip()
        elif line.strip() == "" and a_line != None:
            ids = []
            scfs = []
            start_poss = []
            end_poss = []
            seqs = []
            for s_line in s_lines:
                ary = s_line.split()
                ids.append(ary[1].split(".")[0]) # This assumes that IDs are in the format "species.scaffold"
                scfs.append(ary[1].split(".")[1])
                start_poss.append(int(ary[2]))
                end_poss.append(int(ary[2]) + int(ary[3]))
                seqs.append(ary[6])
            # Use only alignments with sufficient sequences.
            if len(ids) >= min_n_seqs:
                # Make sure that the ids are unique to exclude paralogs.
                if len(set(ids)) == len(ids):
                    # Make sure that the ids are unchanged from the first alignment.
                    if first_set_of_ids == None:
                        first_set_of_ids = ids
                        for idx in first_set_of_ids:
                            vcf_string += "\t" + idx
                        vcf_string += "\n"
                    else:
                        if ids != first_set_of_ids:
                            print("ERROR: The order of IDs seems to vary between alignments!")
                            sys.exit(0)
                    # Make sure that all sequences have the same length.
                    for seq in seqs:
                        if len(seq) != len(seqs[0]):
                            print("ERROR: Sequences differ in length!")
                            sys.exit(1)
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
                            # Identify polymorphic sites in the alignment.
                            for pos in range(len(seqs[0])):
                                non_missing_chars_this_site = []
                                for seq in seqs:
                                    if seq[pos] in ["A", "a", "C", "c", "G", "g", "T", "t"]:
                                        non_missing_chars_this_site.append(seq[pos].upper())
                                # Prepare a vcf line for this polymorphic site.
                                if len(set(non_missing_chars_this_site)) > 1:
                                    all_chars_this_site = []
                                    for seq in seqs:
                                        all_chars_this_site.append(seq[pos].upper())
                                    ref_allele = all_chars_this_site[0] # This assumes that the first sequence can be used as the reference.
                                    if ref_allele in ["A", "C", "G", "T"]:
                                        alt_alleles = []
                                        for char in all_chars_this_site:
                                            if char != ref_allele:
                                                if char in ["A", "C", "G", "T"]:
                                                    if char not in alt_alleles:
                                                        alt_alleles.append(char)
                                        # Prepare a string for the alt alleles (e.g. "A,T")
                                        alt_alleles_as_string = ""
                                        for alt_allele  in alt_alleles:
                                            alt_alleles_as_string += alt_allele + ","
                                        alt_alleles_as_string = alt_alleles_as_string[:-1]
                                        # Get the genotype indices.
                                        gt_indices = []
                                        for char in all_chars_this_site:
                                            if char == ref_allele:
                                                gt_indices.append("0")
                                            elif char in ["A", "C", "G", "T"]:
                                                gt_indices.append(str(alt_alleles.index(char)+1))
                                            else:
                                                gt_indices.append(".")
                                        # Compile the vcf line.
                                        vcf_string += "\t".join([scfs[0], str(start_poss[0]+pos), ".", ref_allele, alt_alleles_as_string, ".", "PASS", ".", "GT"])
                                        for gt_index in gt_indices:
                                            vcf_string += "\t" + gt_index + "|" + gt_index
                                        vcf_string += "\n"
                                        n_variants += 1
            a_line = None
            s_lines = []
        elif line[0:2] == "s\t" and a_line != None:
            s_lines.append(line.strip())

# Write the vcf output file.
with open(vcf_name, "w") as f:
    f.write(vcf_string)
    print("Wrote file " + vcf_name + " with " + str(n_variants) + " variants.")
