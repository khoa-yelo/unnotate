import argparse
import os
from os.path import join
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import esm
from .protein_embedder import embed_proteins
from .faiss_knn import FaissKNN
from .pairwise_alignment_fast import load_fasta_sequences, calculate_sequence_identity
import zipfile

def parse_args():
    parser = argparse.ArgumentParser(description="Unnotate protein sequences")
    parser.add_argument("--fasta_file", type=str, required=True, help="Path to the FASTA file")
    parser.add_argument("--embeddings_db", type=str, required=True, help="Path to the HDF5 embeddings database file") 
    parser.add_argument("--uniprot_db", type=str, required=True, help="Path to the CSV uniprot database file")
    parser.add_argument("--k", type=int, default=20, help="Number of nearest neighbors to search")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to the output directory")
    parser.add_argument("--metric", type=str, default="mean_middle_layer_12", help="Metric to use for the embeddings")
    parser.add_argument("--prefix", type=str, default="unnotated", help="Prefix for output files")
    return parser.parse_args()

def load_embeddings(embeddings_db, metric):
    """Load embeddings from HDF5 file and select metric."""
    with h5py.File(embeddings_db, "r") as f:
        embeddings = {
            "mean": f["mean"][...],
            "max": f["max"][...],
            "mean_middle_layer_12": f["mean_middle_layer_12"][...],
            "max_middle_layer_12": f["max_middle_layer_12"][...]
        }
        data_ids = f["indices"][...]
    
    if metric not in embeddings:
        raise ValueError(f"Invalid metric: {metric}")
        
    return embeddings[metric], data_ids

def find_nearest_neighbors(query_embeddings, reference_embeddings, k, data_ids):
    """Find k nearest neighbors using Faiss."""
    knn = FaissKNN(dim=reference_embeddings.shape[1], metric="cosine", use_gpu=True)
    knn.add(reference_embeddings)
    similarity, indices = knn.search(query_embeddings, k=k)
    
    similarity = similarity.flatten()
    indices = indices.flatten()
    accessions = data_ids[indices].astype(str).flatten()
    
    return similarity, indices, accessions

def save_results(output_dir, df, accessions, similarity, percent_identity_matrix, prefix="unnotated"):
    """Save results to output directory and create zip file."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Save individual files
    csv_path = join(output_dir, f"{prefix}_uniprot.csv")
    acc_path = join(output_dir, f"{prefix}_accessions.npy") 
    dist_path = join(output_dir, f"{prefix}_distances.npy")
    ident_path = join(output_dir, f"{prefix}_sequence_identity.npy")
    
    df.to_csv(csv_path, index=False)
    np.save(acc_path, accessions)
    np.save(dist_path, similarity)
    np.save(ident_path, percent_identity_matrix)
    
    # Create zip file containing all outputs
    zip_path = join(output_dir, f"{prefix}_results.zip")
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        zipf.write(csv_path, arcname=f"{prefix}_uniprot.csv")
        zipf.write(acc_path, arcname=f"{prefix}_accessions.npy")
        zipf.write(dist_path, arcname=f"{prefix}_distances.npy") 
        zipf.write(ident_path, arcname=f"{prefix}_sequence_identity.npy")
    
    print(f"Saved sequence identity matrix to {ident_path}")
    print(f"Created zip archive of all results at {zip_path}")

def calculate_sequence_identities(query_seqs, k, accessions_2d, acc_to_seq):
    """Calculate sequence identity between queries and their matches."""
    N_query = len(query_seqs)
    percent_identity_matrix = np.zeros((N_query, k), dtype=np.float32)
    
    for i, query_seq in enumerate(query_seqs):
        for j in range(k):
            acc = accessions_2d[i, j]
            target_seq = acc_to_seq.get(acc, "")
            percent_identity_matrix[i, j] = calculate_sequence_identity(query_seq, target_seq)
            
    return percent_identity_matrix

def unnotate(fasta_file, embeddings_db, uniprot_db, k=20, metric="mean_middle_layer_12", output_dir=None, prefix="unnotated"):
    # Load embeddings and find reference data
    reference_embeddings, data_ids = load_embeddings(embeddings_db, metric)
    
    # Generate query embeddings
    protein_embeddings = embed_proteins(protein_sequences=[], fasta_file=fasta_file)
    query_embeddings = protein_embeddings[metric]
    
    # Find nearest neighbors
    similarity, indices, accessions = find_nearest_neighbors(
        query_embeddings, reference_embeddings, k, data_ids
    )
    
    # Load and filter UniProt data
    df = pd.read_csv(uniprot_db)
    df = df[df["accession"].isin(accessions)]
    
    if output_dir:
        # Load sequences and prepare for alignment
        print("Performing sequence alignment of queries to their k nearest neighbors...")
        query_seqs = load_fasta_sequences(fasta_file)
        query_seqs = [q[1] for q in query_seqs]
        
        # Build mapping from accession to sequence
        acc_to_seq = {row["accession"]: row["sequence"] 
                     for _, row in df.iterrows() 
                     if pd.notna(row["sequence"])}
        
        # Reshape indices for sequence identity calculation
        indices_2d = indices.reshape(-1, k)
        accessions_2d = data_ids[indices_2d].astype(str)
        
        # Calculate sequence identities
        percent_identity_matrix = calculate_sequence_identities(
            query_seqs, k, accessions_2d, acc_to_seq
        )
        
        # Save all results
        save_results(output_dir, df, accessions, similarity, percent_identity_matrix, prefix)

if __name__ == "__main__":
    args = parse_args()
    unnotate(args.fasta_file, args.embeddings_db, args.uniprot_db, args.k, args.metric, args.output_dir, args.prefix)