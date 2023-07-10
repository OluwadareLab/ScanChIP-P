#######################################
#              ScanChIP-P             #
#         By: Ashley Doerfler         #
#######################################

from sklearn.cluster import DBSCAN
import pandas as pd
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
from knee import KneeLocator

def main():
    # read HP contact matrix
    # hp_matrix_file = f"../examples/hp_matrix.txt" #use parser for this?
    # real_data = read_file(hp_matrix_file, "/t")
    args_length = len(args)
    inputfilename = args[0]
    Resolution = int(args[1])

    if args_length == 2:
        pass
    elif args_length == 3:
        window = args.window
    elif args_length == 4:
        window = int(args[2])
        min_TAD_size = int(args[3])
    elif args_length == 5:
        window = int(args[2])
        min_TAD_size = int(args[3])
        max_TAD_size = int(args[4])
    # Split inputfilename and extract the file name
    tmp = inputfilename.split("[/\\.\\\\ ]")
    if "." in inputfilename:
        name = tmp[-2]
    else:
        name = tmp[-1]
    # Set output paths
    Outputpath = "Output_" + name + "/"
    Clusterpath = "Clusters/"
    ClusterFolder = os.path.join(Outputpath, Clusterpath)
    TADFolder = os.path.join(Outputpath, "TADs/")
    fname = inputfilename
    # create new feature data
    feat = create_new_data(real_data)
    # Create an output folder
    os.makedirs(Outputpath, exist_ok=True)
    filename = "Readme.txt"
    file_path = os.path.join(Outputpath, filename)
    # Initialize log_outputWriter
    log_outputWriter = open(file_path, 'w')
        

#######################################
#          Read HP Matrix File        #
#######################################
def read_file(filename, sep):
    rows = 0
    cols = 0
    try:
        with open(filename, 'r') as file:
            for line in file:
                rowdata = line.rstrip('\n')
                line = rowdata.split(sep)
                rows += 1
                cols = len(line)

        print("Number of row/col =", rows)

        a = [[0.0] * cols for _ in range(rows)]
        lines_counter = 0
        with open(filename, 'r') as file:
            for line in file:
                rowdata = line.rstrip('\n')
                line = rowdata.split(sep)
                for k in range(cols):
                    a[lines_counter][k] = float(line[k])
                lines_counter += 1

        global nRegion
        nRegion = rows

        return a
    except FileNotFoundError as e:
        print("File not found:", filename)
        raise e

    

#######################################
#          Feature Extraction         #
#######################################

# Find zero rows from data
def find_zero_rows(data):
    global zero_rows
    zero_rows = []
    num_dimensions = len(data)
    for dim in range(num_dimensions):
        count = 0
        for value in data[dim]:
            if value == 0:
                count += 1
        if count == len(data[dim]):
            zero_rows.append(dim)
            
# Create a window for feature extraction to exclude a lot of the noise
def create_new_data(mat):
    find_zero_rows(mat)
    print("The number of Zero Rows =", len(zero_rows))
    nrows = nRegion - len(zero_rows)
    ncols = 2 * nRegion
    feature = []
    index = 0
    for diag in range(nRegion):
        if diag in zero_rows:
            continue
        else:
            mat.rowdata(mat, diag, ncols)
            mat.coldata(mat, diag, ncols)
            feature.append(list[:ncols])
            index += 1
            list.clear()
    return feature

#######################################
#          DBScan Clustering          #
#######################################

# Calculate ESP with k-distance
def k_distance(matrix, k):
    n = matrix.shape[0]
    distances = np.zeros(n)

    for i in range(n):
        row = matrix[i]
        row_sorted = np.sort(row)
        distances[i] = row_sorted[k]

    return distances

# Calculate min_points (min tad size / bin resolution)
def min_pnts(min_tad_size, resolution):
    return min_tad_size // resolution

#######################################
#            Generate TADs            #
#######################################


#######################################
#        Evaluate TAD Quality         #
#######################################


#######################################
#            Set Up Parser            #
#######################################
def setup_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', help='input directory containing bedpe files', required=True)
    parser.add_argument('-o', '--outdir', help='output directory', required=True)
    parser.add_argument('-p', '--prefix', help='name of the dataset; should be the same prefix that input bedpe files start with')
    parser.add_argument('-C', '--chip', help='filepath to the ChIP peaks', required=True)
    parser.add_argument('-c', '--chromosomes', help='A comma-separated list of chromosomes', required=True)
    parser.add_argument('-f', '--features', help='path of the file containing genomic features (mappability, gccontent, effective length)', required=True)
    parser.add_argument('-x', '--filter', help='None" or path to the bed file containing regions to be filtered', default='None', required=False)
    parser.add_argument('-b', '--binsize', help='bin size, default = 40000', default=40000, type=int)
    parser.add_argument('-u', '--upperlimit', help='upper limit for distance between bins, default = 2000000', default=2000000, type=int)
    return parser

def parse_arguments(parser):
    args = parser.parse_args()
    print('Reading input from directory:', args.indir)
    print('Output will be written to:', args.outdir)
    print('Genomic features file specified:', args.features)
    chromosomes = ['chr' + chrom for chrom in args.chromosomes.split(',')]
    args.chroms = chromosomes
    print('Analyzing')
    print('   ', args.chroms)
    if args.filter == 'None':
        print('No filter region specified')
    else:
        print('Filtered bins will be read from:', args.filter)
    print('Binsize:', args.binsize)
    print('Prefix:', args.prefix)
    print('Upperlimit:', args.upperlimit)
    return args