#######################################
#              ScanChIP-P             #
#         By: Ashley Doerfler         #
#######################################

from sklearn.cluster import DBSCAN
from sklearn.cluster import HDBSCAN
from sklearn.neighbors import NearestNeighbors
import pandas as pd
import numpy as np
from numpy import savetxt
import seaborn as sns
from sklearn.metrics import silhouette_score
from jqmcvi import base
import argparse
from sklearn.cluster import KMeans
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
import re

def main():
    """
    Read inputs and creates output TADs
    """
    
    parser = setup_parser()
    args = parse_arguments(parser)
    
    # min_points = min_pnts(args.minsize, args.binsize)
    
    print("Read data...")
    data = pd.read_csv(args.input, header=None, sep='\t') # remove headers as needed
    
    # data = submatrix(data, 0, 0, 100) # For testing purposes on simulated data from CASPIAN
    
    data = data.to_numpy()
    
    # np.set_printoptions(precision=3, suppress=True)
    # with open('100x100.txt', 'w') as outfile: # uncomment if using simulated matrix
    #     np.savetxt(outfile, data)
    
    new_matrix = data

    for i in range(0, 30):
        for j in range(i, 30):
            new_matrix[i, j] = data[j,i]
        
    print(new_matrix)
    
    savetxt('sym_30x30.txt', new_matrix, delimiter=' ')
    
    window_length = len(new_matrix) // args.windowproportion
    
    print("Create features...")
    features = create_feature_data(new_matrix, window_length, args.windowproportion)
    
    # Save features to a file
    with open('features2.txt', 'w') as outfile:
        np.savetxt(outfile, features)
    
    print("Find best eps from elbow in range of 10 and cluster...")
    # min_points = min_pnts(args.minsize, args.binsize)
    best_clusters = cluster_ranges(eps=29, features=features)
    
    tads = generate_tad(bin_size=args.binsize, tad_size=args.minsize, clusters=best_clusters)
    
    print(tads)
    tadFile = 'TADfile.txt'
    
    tadQuality(tadFile=tadFile, data=data)
    
    create_heatmap(hic=data, tadfile=tadFile)

def submatrix( matrix, startRow, startCol, size):
    x = np.array(matrix)
    return x[startRow:startRow+size,startCol:startCol+size]

#######################################
#          Feature Extraction         #
#######################################
            
def create_feature_data(matrix, window_length, windowproportion):
    """
    Create a window for feature extraction to exclude a lot of the noise

    Args:
        matrix (np matrix): Contact Matrix
        window (float): Portion of the data to focus on

    Returns:
        array: Features, N_features = (2 * window_length) * len(matrix)
    """
    # Cross Features
    features = []
    window_length = len(matrix) - windowproportion
    
    for diag in range(0, len(matrix)):
        feat = []
        # Finds end position of window
        window = window_length // 2
        
    # Collect row and column features
        # Top-left corner
        if (diag - window) < 0:
            for j in range(0, window_length):
                # row
                feat.append(matrix[diag, j])
            for i in range(0, window_length):
                # col
                feat.append(matrix[i, diag])
        # Bottom-right corner
        elif (diag + window) > len(matrix):
            for j in range(len(matrix) - window_length, len(matrix)):
                feat.append(matrix[diag, j])
            for i in range(len(matrix) - window_length, len(matrix)):
                feat.append(matrix[i, diag])
        # Full cross, and window length is an even nubmer
        elif window_length % 2 == 0:   
            for j in range(diag - window, diag + window): 
                # rows
                feat.append(matrix[diag, j])
            for i in range(diag - window, diag + window):
                # col
                feat.append(matrix[i, diag])
        # Window is not even and in bottom-right corner
        elif diag > (len(matrix) - window_length):
            for j in range((diag - window) - 1, diag + window): 
                # rows
                feat.append(matrix[diag, j])
            for i in range((diag - window) - 1, diag + window):
                # col
                feat.append(matrix[i, diag])
        # Window is not even and in top-left corner
        else:
            for j in range(diag - window, (diag + window) + 1): 
                # rows
                feat.append(matrix[diag, j])
            for i in range(diag - window, (diag + window) + 1):
                # col
                feat.append(matrix[i, diag])
                
        # Add the feature to the set of features            
        features.append(feat)
        
    # print(features)
    features = np.array(features, dtype=float)
    print(features)
    return features

    # L-Shaped Features
    # features = []
    # for diag in range(0, len(matrix)):

    #      # Finds end position of window
    #     window = window_length + diag

    #  # Collect row and column features

    #      # If window is out of bounds
    #     if (window > len(matrix)):

    #          # row
    #         for j in range(0, window_length):
    #              features.append(matrix[diag, diag - j])
    #         # col
    #         for i in range(0, window_length):
    #             features.append(matrix[diag - i, diag])
    #     else:   
    #         for j in range(diag, window): 
    #             # rows
    #             features.append(matrix[diag, j])
    #         for i in range(diag, window):
    #              # col
    #              features.append(matrix[i, diag])

    # return features

    # Square Features
    # features = []
    # for diag in range(0, len(matrix)):
    #     feat = []
    #     # Finds end position of window
    #     window = window_length + diag
        
    # # Collect row and column features

    #     # If window is out of bounds
    #     if (window > len(matrix)):
    #         for i in range(0, window_length): 
    #             # rows
    #             for j in range(0, window_length):
    #                 # col
    #                 feat.append(matrix[diag - i, diag - j])
    #     else:   
    #         for i in range(diag, window): 
    #             # rows
    #             for j in range(diag, window):
    #                 # col
    #                 feat.append(matrix[i, j])
    #     # Add the feature to the set of features    
    #     features.append(feat)
    
    # features = np.array(features, dtype=float)
    
    # return features
    
    # ClusterTAD Features
    # features = []
    # n = 0

    # for diag in range(0, len(matrix)):
    #     feat = []
    #     # Finds end position of window
    #     window = window_length + diag
        
    # # Collect row and column features 
    #     # if diag > len(matrix) / 2:
    #     #     for j in range(n, len(matrix)): 
    #     #     # rows
    #     #         feat.append(matrix[diag, j])
    #     #     for i in range(n, len(matrix)):
    #     #     #      # col
    #     #         feat.append(matrix[i, diag])
    #     # else:
    #     for j in range(0, len(matrix) - n): 
    #             # rows
    #         feat.append(matrix[diag, j])
    #     for i in range(0, len(matrix) - n):
    #         #      # col
    #         feat.append(matrix[i, diag])
                
    #     # Add the feature to the set of features            
    #     features.append(feat)
        
    # # print(features)
    # features = np.array(features, dtype=float)
    
    # print(features)
    # return features

#######################################
#          DBScan Clustering          #
#######################################

def knn(min_pnts, features):
    """Calculate K-NearestNeighbor

    Args:
        min_pnts (_type_): _description_
        features (_type_): _description_

    Returns:
        _type_: _description_
    """
    neighbors = NearestNeighbors(n_neighbors=min_pnts)
    neighbors_fit = neighbors.fit(features)
    distances, indices = neighbors_fit.kneighbors(features)
    distances = np.sort(distances, axis=0)
    distances = distances[:,1]
    plt.plot(distances) 
    plt.savefig('KNN2.png')
    
    return distances

def find_elbow_point(k_distances):
    """
    Find the elbow point in the k-distances curve using the elbow method.

    Args:
        k_distances (array): The k-distances obtained from the k_distance function.

    Returns:
        int: The index of the elbow point in the k-distances array.
    """
    distortions = []
    for i in range(1, len(k_distances)):
        # Calculate the squared distance between each k-distance and its previous one
        distortion = (k_distances[i] - k_distances[i - 1]) ** 2
        distortions.append(distortion)

    # Find the index of the elbow point where the slope starts to decrease significantly
    elbow_point_index = distortions.index(max(distortions))
    print(elbow_point_index)
    return elbow_point_index + 1  # Add 1 to get the actual k value for the elbow point

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

def cluster_ranges(eps, features):
    """Finds clusters for eps ranging from the estimate eps -5 to +5

    Args:
        eps (int): DBSCAN radius
        min_points (int): DBSCAN minimum neighbors
        features (2D array): Dataset for Clustering
    """
    
    best_score_s = 0
    best_file = ""
    best_cluster = []
    for e in range(1, 10):
        
        # Cluster TADs
        # clusters = DBSCAN(eps=e, min_samples=3).fit(features)
        # clusters = HDBSCAN().fit(features)
        clusters = KMeans(n_clusters=e, random_state=0, n_init="auto").fit(features)
        
        # Count the number of different labels
        labels = clusters.labels_
        print(labels)
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            
        if n_clusters_ <= 1:
            continue
            
        filename = f"cluster_eps{e}_minpoints_2.txt"
            
        # Silhouette Score
        s_s = silhouette_score(features, labels)
        if s_s > best_score_s:
            best_score_s = s_s
            best_file = filename
            best_cluster = clusters.labels_
        print("Silhouette Score(n=3): ", s_s)
                
    print("The best file was ", best_file, "with a Silhouette Score of ", best_score_s)
    print(best_cluster)
    
    return best_cluster

#######################################
#            Generate TADs            #
#######################################

def generate_tad(bin_size, tad_size, clusters):
    """
    Generate TADs based on how many bins are required for the minimum size of a TAD
    """
    
    count = 0
    start = None
    tad_count = 0
    previous = None
    tads = []
    
    for label in clusters:
        if count >= len(clusters) - 1:
            if tad_count > min_pnts(tad_size, bin_size):
                tads.append([start, count])
        elif previous is None:
            previous = label
            start = count
        elif previous != label:
            if tad_count > min_pnts(tad_size, bin_size):
                tads.append([start, count - 1])
            start = count
            tad_count = 1
        else:
            tad_count += 1
        previous = label
        count +=1
        
    with open('TADfile.txt', 'w') as outfile:
        np.savetxt(outfile, tads)
            
    return tads
    

#######################################
#        Evaluate TAD Quality         #
#        Adapted from CASPIAN         #
#######################################

def tadQuality(tadFile, data):
    """TAD quality"""
    n = len(data)
    tad = np.loadtxt(tadFile, dtype=int)
    num_tads = len(tad)

    intra = 0
    intra_num = 0
    for i in range(num_tads):
        for j in range(int(tad[i, 0]), int(tad[i, 1] + 1)):
            for k in range(int(tad[i, 0]), int(tad[i, 1] + 1)):
                intra += data[j, k]
                intra_num += 1

    if intra_num != 0:
        intra = intra / intra_num
        print("intra TAD: %0.3f" % intra)
    else:
        intra = 0

    inter = 0
    inter_num = 0
    for i in range(num_tads - 1):
        for j in range(int(tad[i, 0]), int(tad[i, 1] + 1)):
            for k in range(int(tad[i + 1, 0]), int(tad[i + 1, 1] + 1)):
                inter += data[j, k]
                inter_num += 1

    if inter_num != 0:
        inter = inter / inter_num
        print("inter TAD: %0.3f" % inter)
    else:
        inter = 0

    print("quality: %0.3f" % (intra - inter))
    quality = (intra - inter) / num_tads
    print(quality)
    
    return quality

#######################################
#       Create Heatmap of TADs        #
#######################################
def readTAD(tadfile):
    f = open(tadfile)
    line=f.readline()
    start=[]
    end=[]
    while line:
        line = line.split()
        try:
            start1 = int(float(line[0]))
        except:
            start1 = 0
        end1 = int(float(line[1]))
        start.append(start1)
        end.append(end1)
        line = f.readline()
    f.close()
    return start, end
# def readTAD(tadfile):
#     f = open(tadfile)
#     line = f.readline()
#     start = []
#     end = []
#     while line:
#         line = line.split()
#         try:
#             start1 = int(float(line[0]))
#         except:
#             start1 = 0
#         end1 = int(float(line[1]))
#         start.append(start1)
#         end.append(end1)
#         line = f.readline()
#     f.close()

#     # Split the start and end lists into separate lists
#     start = start[::2]  # Start from the beginning and take every second element
#     end = end[1::2]     # Start from the second element and take every second element

#     return start, end

def create_heatmap(hic, tadfile):

    # Read TAD boundaries using the readTAD function
    start, end = readTAD(tadfile)
    print("Length of TADs", len(start))
    lentad = len(start)
    tad_label = start + end
    tad_label.sort()

    # Plot the Hi-C matrix heatmap and highlight TAD boundaries
    palette = sns.color_palette("bright", 10)
    plt.figure(figsize=(10.5, 10))
    start1 = 0
    end1 = 7
    try:
        sns.heatmap(data=hic, robust=True, cmap="OrRd")
    except UserWarning:
        sns.heatmap(data=hic[start1:end1, start1:end1], robust=True, cmap="OrRd", vmin=0, vmax=1)

    for i in range(0, lentad):
        if start1 < start[i] < end1 and start1 < end[i]:
            plt.hlines(y=start[i] - start1, xmin=start[i] - start1, xmax=end[i] - start1)
            plt.vlines(x=end[i] - start1, ymin=start[i] - start1, ymax=end[i] - start1)
    plt.title('TAD boundary')
    plt.savefig('../example/{}_{}_{}.pdf'.format(start1, end1, ''), format='pdf', bbox_inches='tight')
    plt.show()

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








