# m_matschiner Sat Dec 4 00:18:33 CET 2021

# Load libraries.
import numpy as np
import pandas as pd
import sys

# Get the command-line arguments.
in_file_name = sys.argv[1]
out_file_name = sys.argv[2]

# Prepare an empty square matrix.
num_lines_in_file = sum(1 for line in open(in_file_name))
distances = []
sample_names = []

# Read the input file.
with open(in_file_name) as f:
    next(f) # skip sample count line
    for line in f:
        elements = line.strip().split('\t')
        sample_names.append(elements[0])
        row = [float(e) for e in elements[1:]]
        row.extend([0.0] * (num_lines_in_file-1-len(row)))
        distances.append(row)
    np_array = np.asarray(distances)
    index_upper = np.triu_indices(num_lines_in_file-1)
    np_array[index_upper] = np_array.T[index_upper]
    matrix = pd.DataFrame(np_array, columns=sample_names, index=sample_names).to_csv(sep='\t',header=False)

# Write the output file.
with open(out_file_name, "w") as out_file:
    out_file.write(str(len(sample_names)) + "\n")
    out_file.write(matrix)
