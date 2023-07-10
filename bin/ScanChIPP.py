#######################################
#              ScanChIP-P             #
#         By: Ashley Doerfler         #
#            Date: 07-2023            #
#######################################

from sklearn.cluster import DBSCAN
import pandas as pd
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
from knee import KneeLocator

def main():
    """
    Read inputs and creates output TADs
    """
    # read HP contact matrix
    real_data = read_file(args.input, "/t")
    Resolution = args.binsize
    window = int(args.window)
    min_TAD_size = args.min_tad_size
    max_TAD_size = args.max_tad_size
    # Split inputfilename and extract the file name
    tmp = args.input.split("[/\\.\\\\ ]")
    if "." in args.input:
        name = tmp[-2]
    else:
        name = tmp[-1]
    # Set output paths
    Outputpath = "Output_" + name + "/"
    Clusterpath = "Clusters/"
    ClusterFolder = os.path.join(Outputpath, Clusterpath)
    TADFolder = os.path.join(Outputpath, "TADs/")
    fname = arg.input
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
    """
    Read the contact matrix file

    Args:
        filename (String): Name of the file that needs to be read
        sep (String): How to seperate

    Raises:
        e: File not found

    Returns:
        Matrix: Parsed contact matrix
    """
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

        contact_matrix = [[0.0] * cols for _ in range(rows)]
        lines_counter = 0
        with open(filename, 'r') as file:
            for line in file:
                rowdata = line.rstrip('\n')
                line = rowdata.split(sep)
                for k in range(cols):
                    contact_matrix[lines_counter][k] = float(line[k])
                lines_counter += 1

        global nRegion
        nRegion = rows

        return contact_matrix
    except FileNotFoundError as e:
        print("File not found:", filename)
        raise e

    

#######################################
#          Feature Extraction         #
#######################################

def find_zero_rows(matrix):
    """
    Find zero rows from data

    Args:
        matrix (np matrix): Contact Matrix
    """
    global zero_rows
    zero_rows = []
    num_dimensions = len(matrix)
    for dim in range(num_dimensions):
        count = 0
        for value in matrix[dim]:
            if value == 0:
                count += 1
        if count == len(matrix[dim]):
            zero_rows.append(dim)
            
def create_new_data(matrix):
    """
    Create a window for feature extraction to exclude a lot of the noise

    Args:
        matrix (np matrix): Contact Matrix

    Returns:
        _type_: _description_
    """
    find_zero_rows(matrix)
    print("The number of Zero Rows =", len(zero_rows))
    nrows = nRegion - len(zero_rows)
    ncols = 2 * nRegion
    feature = []
    index = 0
    for diag in range(nRegion):
        if diag in zero_rows:
            continue
        else:
            matrix.rowdata(matrix, diag, ncols)
            matrix.coldata(matrix, diag, ncols)
            feature.append(list[:ncols])
            index += 1
            list.clear()
    return feature

#######################################
#          DBScan Clustering          #
#######################################

def k_distance(matrix, k):
    """
    Calculate ESP with k-distance

    Args:
        matrix (_type_): _description_
        k (_type_): _description_

    Returns:
        _type_: _description_
    """
    n = matrix.shape[0]
    distances = np.zeros(n)

    for i in range(0, n - 1):
        row = matrix[i]
        row_sorted = np.sort(row)
        distances[i] = row_sorted[k]

    return distances

def min_pnts(min_tad_size, resolution):
    """
    Calculate min_points ~or~ min_bins (min tad size / bin resolution)

    Args:
        min_tad_size (_type_): Minimum size of a TAD default=120k
        resolution (_type_): Binsize default=40kb

    Returns:
        int: returns the amout of bins required to be considered a TAD
    """
    return min_tad_size // resolution

#######################################
#            Generate TADs            #
#######################################

def generate_tad():
    """
    Generate TADs based on how many bins are required for the minimum size of a TAD
    """
    return

#######################################
#        Evaluate TAD Quality         #
#######################################

def measure_of_concordance(tad_1, tad_2):
    """
    Computes the Measure of Concordance (MoC) between two sets of TAD regions.

    Parameters:
    tad_1 (list): The first set of TAD regions.
    tad_2 (list): The second set of TAD regions.

    Returns:
    float: The Measure of Concordance (MoC) between A and B.
    """
    n_tad_1 = len(tad_1)
    n_tad_2 = len(tad_2)

    sum_term = 0

    for i in range(n_tad_1):
        for j in range(n_tad_2):
            common_bins = len(set(tad_1[i]) & set(B[j]))
            sum_term += (common_bins ** 2) / (len(tad_1[i]) * len(tad_2[j]) - 1)

    moc = (1 / (n_tad_1 * n_tad_2 - 1)) * sum_term
    return moc


def modified_jaccard_index(tad_1, tad_2):
    """
    Computes the modified Jaccard's index between two sets of TAD boundaries.

    Parameters:
    tad_1 (set): The first set of TAD boundaries.
    tad_2 (set): The second set of TAD boundaries.

    Returns:
    float: The modified Jaccard's index between A and B.
    """
    intersecting_set = set()
    double_counted_boundaries = set()

    for boundary in tad_1:
        for offset in range(-1, 2):
            shifted_boundary = boundary + offset
            if shifted_boundary in tad_2:
                if shifted_boundary not in double_counted_boundaries:
                    intersecting_set.add(shifted_boundary)
                    double_counted_boundaries.add(shifted_boundary)

    union_size = len(tad_1) + len(tad_2) - len(intersecting_set)
    jaccard_index = len(intersecting_set) / union_size if union_size != 0 else 0

    return jaccard_index


def length_quality():
    """
    Find the quality of the TAD predictions based on length of the TADs
    """
    return

def amount_identified_quality():
    """
    Find the quality of the TAD predictions based on the number of identified TAD's
    """
    return

def interaction_quality():
    """
    Intra- and inter- Cluster (TAD) similarity
    """
    return

#######################################
#            Set Up Parser            #
#######################################
def setup_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input contact matrix files', required=True)
    parser.add_argument('-w', '--window', help='window', required=True)
    parser.add_argument('-m', '--minsize', help='minimum size of a TAD', default=120000, type=int)
    parser.add_argument('-b', '--binsize', help='bin size, default = 40000', default=40000, type=int)
    return parser

def parse_arguments(parser):
    args = parser.parse_args()
    print('Contact matrix specified:', args.input)
    print('Creating window:', args.window)
    print('Binsize:', args.binsize)
    print('Minimum TAD size', args.minsize)
    print('Maximum TAD size', args.maxsize)
    return args