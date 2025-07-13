import pickle
import h5py
import numpy as np

path = "/home/jovyan/workspace/data/unnot/parsed_uniprot_swiss_data_esm_input.pkl"
with open(path, 'rb') as f:
    data = pickle.load(f)

metrics = ['mean', 'max', 'mean_middle_layer_12', 'max_middle_layer_12']
indices = list(data.keys())
num_items = len(indices)

# Infer embedding dimension from the first entry
first_key = indices[0]
dims = {metric: data[first_key][metric].shape[0] for metric in metrics}

with h5py.File('embeddings.h5', 'w') as f:
    # Create datasets for each metric
    dsets = {}
    for metric in metrics:
        dsets[metric] = f.create_dataset(
            metric,
            shape=(num_items, dims[metric]),
            dtype='float32'
        )
    # Create dataset for indices (as fixed-length ASCII strings)
    max_index_len = max(len(idx) for idx in indices)
    dset_indices = f.create_dataset(
        'indices',
        shape=(num_items,),
        dtype=f'S{max_index_len}'
    )

    # Fill datasets
    for i, idx in enumerate(indices):
        for metric in metrics:
            dsets[metric][i] = data[idx][metric]
        dset_indices[i] = np.string_(idx)

print("HDF5 file 'embeddings.h5' created with datasets for each metric and indices.")