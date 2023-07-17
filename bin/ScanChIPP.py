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

def main():
    """
    Read inputs and creates output TADs
    """
    
    # Setup parser
    parser = setup_parser()
    args = parse_arguments(parser)
    
    # read HP contact matrix
    print("Read data")
    data = pd.read_csv(args.input, header=None, sep='\t')
    data = data.to_numpy()
    # print(data[1, :])
    
    # Create feature data from contact matrix with specified window size
    print("Create features")
    features = create_feature_data(data, args.windowproportion, args.binsize)
    # what do the features look like
    
    # print("Make clusters")
    # clusters = DBSCAN(eps=k_distance(features), min_samples=min_pnts(args.minsize, args.binsize)).fit(features)
    # print(clusters.labels_) # just returns DBSCAN(eps=1, min_samples=3) ????

#######################################
#          Feature Extraction         #
#######################################
            
def create_feature_data(matrix, window_proportion):
    """
    Create a window for feature extraction to exclude a lot of the noise

    Args:
        matrix (np matrix): Contact Matrix
        window (float): Portion of the data to focus on

    Returns:
        array: Features, N_features = (2 * window_length) * len(matrix)
    """
    features = []

    for diag in range(0, len(matrix)):

        # Finds end position of window
        window = (len(matrix) // window_proportion) + diag
        
        # Collect row and column features

        # If window is out of bounds
        if (window > len(matrix)):
            window_length = len(matrix) // window_proportion
            # row
            for j in range(0, window_length):
                features.append(matrix[diag, diag - j])
            # col
            for i in range(0, window_length):
                features.append(matrix[diag - i, diag])
        else:   
            for j in range(diag, window): 
                # rows
                features.append(matrix[diag, j])
            for i in range(diag, window):
                # col
                features.append(matrix[i, diag])

    print(features, "length: ", len(features)) 
            
    return features

#######################################
#          DBScan Clustering          #
#######################################

def k_distance(matrix):
    """
    Calculate ESP with k-distance

    Args:
        matrix (_type_): _description_
        k (_type_): _description_

    Returns:
        _type_: _description_
    """
    n = len(matrix)
    distances = np.zeros((n, n-1))

    for i in range(n):
        row = matrix[i]
        row_sorted = np.sort(row)
        for k in range(1, n):
            distances[i, k-1] = row_sorted[k]

    elbow = find_elbow(distances)
    return elbow

def find_elbow(distances):
    """
    Find the elbow point in a curve.

    Args:
        distances (ndarray): Array of distances.

    Returns:
        int: Index of the elbow point.
    """
    n = len(distances)
    x = np.arange(1, n + 1)  # Indices of the distances

    # Compute the cumulative sum of squared differences
    cum_sum_sq_diff = np.cumsum((distances - distances.mean())**2)

    # Compute the total sum of squared differences
    total_sq_diff = cum_sum_sq_diff[-1]
    
    # Check if total_sq_diff is zero
    if total_sq_diff == 0:
        explained_variance = np.zeros_like(cum_sum_sq_diff)
    else:
        # Replace zeros in total_sq_diff with a small nonzero value
        total_sq_diff = np.where(total_sq_diff == 0, 1e-10, total_sq_diff)
        # Compute the ratio of explained variance
        explained_variance = cum_sum_sq_diff / total_sq_diff
    # Compute the ratio of explained variance
    explained_variance = cum_sum_sq_diff / total_sq_diff

    # Find the elbow point
    elbow_index = np.argmax(explained_variance >= 0.9) + 1

    # Plot the explained variance curve
    plt.plot(x, explained_variance, marker='o')
    plt.xlabel('Number of clusters')
    plt.ylabel('Explained variance')
    plt.title('Explained Variance Curve')
    plt.show()

    return elbow_index

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
    parser.add_argument('-i', '--input', help='input contact matrix file', required=True)
    parser.add_argument('-w', '--windowproportion', help='window proportion', type=int)
    parser.add_argument('-m', '--minsize', help='minimum size of a TAD', default=120000, type=int)
    parser.add_argument('-b', '--binsize', help='bin size, default = 40000', default=40000, type=int)
    return parser

def parse_arguments(parser):
    args = parser.parse_args()
    print('Contact matrix specified:', args.input)
    print('Creating window:', 1 / args.windowproportion)
    print('Binsize:', args.binsize)
    print('Minimum TAD size', args.minsize)
    return args

main()

