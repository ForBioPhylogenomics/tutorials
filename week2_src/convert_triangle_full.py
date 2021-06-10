#function from https://github.com/marbl/Mash/issues/9
import numpy as np
import pandas as pd
import sys


def lower_triangle_to_full_natrix(filename):
    num_lines_in_file = sum(1 for line in open(filename))
    distances = []
    sample_names = []

    with open(filename) as f:
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
        print(len(sample_names))
        #print(np_array)
        #return pd.DataFrame(np_array, columns=sample_names, index=sample_names).to_csv('output.tsv', sep="\t")
        for index, i in enumerate(np_array):
            #print(i)
            print(sample_names[index], end='\t')
            print("%s" % '\t'.join([str(x) for x in i.tolist()]))

if __name__ == '__main__':
	lower_triangle_to_full_natrix(sys.argv[1])
