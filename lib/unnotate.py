import argparse
import os
from os.path import join
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import esm
import torch
from .protein_embedder import embed_proteins
from .faiss_knn import FaissKNN
from .pairwise_alignment_fast import load_fasta_sequences, calculate_sequence_identity
import zipfile
import glob
from tqdm import tqdm
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
use_gpu = device.type == "cuda"

def parse_args():
    parser = argparse.ArgumentParser(description="Unnotate protein sequences")
    parser.add_argument("--fasta-file", type=str, required=True, help="Path to the FASTA file")
    parser.add_argument("--database-dir", type=str, required=True, help="Directory containing .h5 (embeddings) and .csv (UniProt) files")
    parser.add_argument("--k", type=int, default=20, help="Number of nearest neighbors to search")
    parser.add_argument("--output-dir", type=str, required=True, help="Path to the output directory")
    parser.add_argument("--metric", type=str, default="mean_middle_layer_12", help="Metric to use for the embeddings")
    parser.add_argument("--prefix", type=str, default="unnotated", help="Prefix for output files")
    parser.add_argument("--faiss-metric", type=str, default="cosine", choices=["cosine", "euclidean"], help="FAISS metric to use (cosine or euclidean)")
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

def find_nearest_neighbors(query_embeddings, reference_embeddings, k, data_ids, database_dir, metric, faiss_metric="cosine"):
    """Find k nearest neighbors using Faiss."""
    # Check if FAISS index already exists
    index_path = os.path.join(database_dir, f"{metric}_{faiss_metric}.faiss.index")
    
    if os.path.exists(index_path):
        logger.info(f"Loading existing FAISS index from {index_path}")
        knn = FaissKNN.load(index_path, metric=faiss_metric, use_gpu=use_gpu)
    else:
        logger.info("Creating new FAISS index")
        knn = FaissKNN(dim=reference_embeddings.shape[1], metric=faiss_metric, use_gpu=use_gpu)
        knn.add(reference_embeddings)
    
    similarity, indices = knn.search(query_embeddings, k=k)
    
    similarity = similarity
    indices = indices
    accessions = data_ids[indices].astype(str)
    
    return similarity, indices, accessions

def save_results(output_dir, df, accessions, similarity, percent_identity_matrix, prefix="unnotated"):
    """Save results to output directory and create zip file. All arrays are saved in a single .npz file."""
    os.makedirs(output_dir, exist_ok=True)

    # Save individual files
    csv_path = join(output_dir, f"{prefix}_uniprot.csv")
    npz_path = join(output_dir, f"{prefix}_results.npz")

    df.to_csv(csv_path, index=False)
    np.savez(npz_path, accession=accessions, cosine_similarity=similarity, sequence_identity=percent_identity_matrix)

    # Create zip file containing all outputs
    zip_path = join(output_dir, f"{prefix}_streamlit.zip")
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        zipf.write(csv_path, arcname=f"{prefix}_uniprot.csv")
        zipf.write(npz_path, arcname=f"{prefix}_results.npz")

    logger.info(f"Saved all results in {npz_path}")
    logger.info(f"Created zip archive of all results at {zip_path}")

def calculate_sequence_identities(query_seqs, k, accessions_2d, acc_to_seq):
    """Calculate sequence identity between queries and their matches."""
    N_query = len(query_seqs)
    percent_identity_matrix = np.zeros((N_query, k), dtype=np.float32)
    
    for i, query_seq in tqdm(enumerate(query_seqs), total=N_query, desc="Calculating sequence identities"):
        for j in range(k):
            acc = accessions_2d[i, j]
            target_seq = acc_to_seq.get(acc, "")
            percent_identity_matrix[i, j] = calculate_sequence_identity(query_seq, target_seq)
            
    return percent_identity_matrix

def unnotate(fasta_file, database_dir, k=20, metric="mean_middle_layer_12", output_dir=None, prefix="unnotated", faiss_metric="cosine"):
    """
    Annotate protein sequences using a database directory containing embeddings and UniProt CSV.
    Args:
        fasta_file: Path to query FASTA file
        database_dir: Directory containing .h5 (embeddings) and .csv (UniProt) files
        k: Number of nearest neighbors
        metric: Embedding metric
        output_dir: Output directory
        prefix: Output file prefix
        faiss_metric: FAISS metric to use ("cosine" or "euclidean")
    """
    # Find files in database_dir
    h5_files = glob.glob(os.path.join(database_dir, "*.h5"))
    csv_files = glob.glob(os.path.join(database_dir, "*.csv"))
    if not h5_files:
        raise FileNotFoundError(f"No .h5 embeddings file found in {database_dir}")
    if not csv_files:
        raise FileNotFoundError(f"No .csv UniProt file found in {database_dir}")
    embeddings_db = h5_files[0]
    uniprot_db = csv_files[0]
    
    # Load embeddings and find reference data
    reference_embeddings, data_ids = load_embeddings(embeddings_db, metric)
    
    # Generate query embeddings
    protein_embeddings = embed_proteins(protein_sequences=[], fasta_file=fasta_file)
    query_embeddings = protein_embeddings[metric]
    
    # Find nearest neighbors
    similarity, indices, accessions = find_nearest_neighbors(
        query_embeddings, reference_embeddings, k, data_ids, database_dir, metric, faiss_metric=faiss_metric
    )
    
    # Load and filter UniProt data
    df = pd.read_csv(uniprot_db)
    df = df[df["accession"].isin(accessions.flatten().astype(str))]
    
    if output_dir:
        # Load sequences and prepare for alignment
        logger.info("Performing sequence alignment of queries to their k nearest neighbors...")
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
    logger.info(f"Starting Unnotate with args: {args}")
    unnotate(
        args.fasta_file,
        args.database_dir,
        args.k,
        args.metric,
        args.output_dir,
        args.prefix,
        faiss_metric=args.faiss_metric
    )