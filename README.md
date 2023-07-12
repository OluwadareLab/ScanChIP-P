# ScanChIP-P
REU 2023 project

ScanChIP-P is a proposed method to identify TAD boundaries from HiChIP and PLAC-seq, protein-centric chromatin conformation methods.

ScanChIP-P implements DBSCAN, a cluster machine learning algorithm to identify clusters of interactions as TADs.

# Input Data
Normalized data was collected from HPTAD and converted into a contact matrix in .hic format
HPTAD Repository: https://github.com/yunliUNC/HPTAD/tree/main

example/hp_matrix.txt

# How to Run
-i Input Matrix File \
-w Window Size (what portion of the data) \
-m Minimum TAD Size \
-b Bin Size 

python3 ScanChIPP.py \
-i  ../example/hp_matrix.txt\
-w .5 \
-m 120000 \
-b 40000

python3 ScanChIPP.py -i ../examples/hp_matrix.txt -w .5 -m 120000 -b 40000
