# m_matschiner Wed May 19 18:38:31 CEST 2021

# This script reads a VCF file with information for a single chromsome,
# and generates a set of alignments in Phylip format for regions
# sampled evenly across the chromosome.

class Window(object):

	def __init__(self, ids, align_length, align_poss, align_refs, align_alts, align_gts):
		DNA_alphabet = ['A','C','G','T']
		genotype_warning_printed = False
		# Prepare an alignment without variation from randomly chosen nucleotides.
		self.ids = ids
		self.seqs = []
		self.n_variable = len(align_poss)
		seq = []
		for _ in range(align_length):
			seq.append(DNA_alphabet[randint(0, 3)])
		for _ in range(len(self.ids)):
			self.seqs.append(copy.deepcopy(seq))
		# Add the variation from the VCF to the alignment.
		for x in range(len(align_poss)):
			pos_in_align = align_poss[x]-1 - align_starts[align_count]
			if pos_in_align < 0 or pos_in_align > align_length - 1:
				print("ERROR: Unexpected alignment position: " + str(pos_in_align) + "!")
				sys.exit(1)
			# Combine ref and alt allleles into a single array.
			ref_alts = [align_refs[x]]
			for alt in align_alts[x].split(","):
				ref_alts.append(alt)
			# Get the array of genotypes as indices.
			gts_as_indices = []
			for align_gt in align_gts[x]:
				if "/" in align_gt:
					if not genotype_warning_printed:
						print("WARNING: Genotypes do not seem to be phased but will be assumed to be so!")
						genotype_warning_printed = True
					align_gt = align_gt.replace("/","|")
				if align_gt.split("|")[0] == ".":
					gts_as_indices.append(-1)
				else:
					gts_as_indices.append(int(align_gt.split("|")[0]))
				if align_gt.split("|")[1] == ".":
					gts_as_indices.append(-1)
				else:
					gts_as_indices.append(int(align_gt.split("|")[1]))
			# Convert genotypes to alleles using the array of ref and alt alleles.
			gts_as_alleles = []
			for gts_as_index in gts_as_indices:
				if gts_as_index == -1:
					gts_as_alleles.append("N")
				else:
					gts_as_alleles.append(ref_alts[gts_as_index])
			# Insert the genotypes into the alignment.
			for y in range(len(self.ids)):
				self.seqs[y][pos_in_align] = gts_as_alleles[y]

	def write_align(self, f_name):
		# Prepare the alignment string in Phylip format.
		outstring = str(len(self.seqs)) + " " + str(len(self.seqs[0])) + "\n"
		for y in range(len(self.seqs)):
			outstring += self.ids[y] + "  " + "".join(self.seqs[y]) + "\n"
		# Write the output file.
		with open(f_name, "w") as f:
			f.write(outstring)
		print("Wrote file " + f_name + " with " + str(self.n_variable) + " variable sites.")

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse
import textwrap
from random import randint
import copy
from pathlib import Path

# Set up the argument parser.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    ------------------------------------------------------------
      This script reads a VCF file with information for a single
      chromsome, and generates a set of alignments in Phylip
      format for regions sampled evenly across the chromosome.

      Run e.g. with
      python3 make_alignments_from_vcf.py input.vcf output_prefix -n 1000 -l 1000 -c 5000000
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
    dest='n_align',
    help='number of alignments (default: 1)'
        )
parser.add_argument(
    '-l',
    nargs=1,
    type=int,
    default=[1],
    dest='align_length',
    help='Length of alignments in bp (default: 1)'
        )
parser.add_argument(
    '-c',
    nargs=1,
    type=int,
    default=[1],
    dest='chr_length',
    help='Length of chromosome in bp  (default: 1)'
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
    'vcf',
    nargs=1,
    type=str,
    help='name of the input file in VCF format'
    )
parser.add_argument(
    'prefix',
    nargs=1,
    type=str,
    help='prefix for output file names'
    )

# Parse command-line arguments.
args = parser.parse_args()
vcf_name = args.vcf[0]
prefix = args.prefix[0]
chr_length = args.chr_length[0]
n_align = args.n_align[0]
align_length = args.align_length[0]
path = args.path[0] + "/"
path = path.replace("//","/")
Path(path).mkdir(parents=True, exist_ok=True)

# Calculate start and end positions of all alignment windows.
align_centers = []
align_center = int(chr_length/n_align/2)
align_centers.append(align_center)
for _ in range(n_align-1):
	align_center += int(chr_length/n_align)
	align_centers.append(align_center)
align_starts = []
align_ends = []
for align_center in align_centers:
	align_starts.append(int(align_center - align_length/2 + 1))
	align_ends.append(int(align_center + align_length/2))

# Parse the VCF input file line by line.
pos = 0
chrom = None
align_count = 0
align_poss = []
align_refs = []
align_alts = []
align_gts = []
ids = []
with open(vcf_name) as vcf:
	for line in vcf:
		if line[0:2] == "##":
			continue
		elif line[0:1] == "#":
			for sample_id in line.split()[9:]:
				ids.append(sample_id + "_1")
				ids.append(sample_id + "_2")
		else:
			prev_pos = pos
			pos = int(line.split()[1])
			# Skip if the current position is 0 (this can happen when the VCF is written with msprime, but not otherwise).
			if pos == 0:
				continue
			# Quit if the current position is smaller than the previous one.
			if pos <= prev_pos:
				print("ERROR: Positions are not continuously increasing, " + str(pos) + " followed " + str(prev_pos) + "!")
				sys.exit(1)
			# Quit if the chromosome has changed.
			if chrom == None:
				prev_chrom = line.split()[0]
			else:
				prev_chrom = chrom
			chrom = line.split()[0]
			if chrom != prev_chrom:
				print("ERROR: The VCF file contains multiple chromosomes!")
				sys.exit(1)
			# Check if the position is within an alignment window.
			if align_count < n_align:
				if pos-1 < align_starts[align_count]:
					continue
				elif pos-1 > align_ends[align_count]:
					# Generate a new alignment.
					window = Window(ids, align_length, align_poss, align_refs, align_alts, align_gts)
					# Write the new alignment.
					f_name = path + prefix + "_" + str(align_starts[align_count] + 1) + "_" + str(align_ends[align_count] + 1) + ".phy"
					window.write_align(f_name)

					# Reset the array of alignment lines and increase the alignment count.
					align_poss = []
					align_refs = []
					align_alts = []
					align_gts = []
					align_count += 1
				else:
					line_ary = line.split()
					align_poss.append(int(line_ary[1]))
					align_refs.append(line_ary[3])
					align_alts.append(line_ary[4])
					align_gts.append(line_ary[9:])

# If no variable sites were in the VCF after the end of an alignment window, some alignments may not have been written yet.
# Complete any unfinished alignments and write these to a files.
while align_count < n_align:
	# Generate a new alignment.
	window = Window(ids, align_length, align_poss, align_refs, align_alts, align_gts)
	# Write the new alignment.
	f_name = path + prefix + "_" + str(align_starts[align_count] + 1) + "_" + str(align_ends[align_count] + 1) + ".phy"
	window.write_align(f_name)

	# Reset the array of alignment lines and increase the alignment count.
	align_poss = []
	align_refs = []
	align_alts = []
	align_gts = []
	align_count += 1
