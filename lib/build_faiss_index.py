#!/usr/bin/env python3
"""
Simple script to build FAISS index from embeddings.h5 using existing FaissKNN class
"""

import h5py
import numpy as np
from pathlib import Path
import sys

# Add lib to path
sys.path.insert(0, str(Path(__file__).parent / "lib"))

from faiss_knn import FaissKNN

# Load embeddings from HDF5
with h5py.File('../database/embeddings.h5', 'r') as f:
    dataset_name = "mean_middle_layer_12"
    embeddings = f[dataset_name][:]
    print(f"Loaded {len(embeddings)} embeddings with dimension {embeddings.shape[1]}")

# Build and save FAISS index
knn = FaissKNN(dim=embeddings.shape[1], metric='cosine', use_gpu=True)
knn.add(embeddings)
knn.save('../database/faiss_index_cosine.faiss')
print("FAISS index saved to ../database/faiss_index_cosine.faiss") 


knn = FaissKNN.load('../database/faiss_index_cosine.faiss', metric='cosine', use_gpu=True)  # or use_gpu=False for CPU
query_embeddings = embeddings[0:5]
# Query with your embeddings (replace query_embeddings with your actual embeddings)
distances, indices = knn.search(query_embeddings, k=5)
print(f"Top 5 matches: {indices}")
print(f"Similarities: {distances}")