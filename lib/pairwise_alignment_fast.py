#!/usr/bin/env python3
"""
Fast Pairwise Sequence Alignment Script using parasail-python

This script performs pairwise sequence alignment between protein sequences from a FASTA file
and sequences from UniProt accessions, reporting percent identity using parasail for speed.

Input:
- FASTA file containing N protein sequences
- NumPy array accession_array with shape (N, k) containing UniProt accessions
- CSV file with UniProt sequences (parsed_uniprot_swiss_data.csv)

Output:
- Percent identity matrix for each protein sequence against its k accessions
"""

import numpy as np
import pandas as pd
import parasail
import argparse
import sys
from Bio import SeqIO
from typing import List, Tuple, Dict
import logging
from tqdm import tqdm
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_fasta_sequences(fasta_file: str) -> List[Tuple[str, str]]:
    """Load protein sequences from FASTA file."""
    sequences = []
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append((record.id, str(record.seq)))
        logger.info(f"Loaded {len(sequences)} sequences from {fasta_file}")
        return sequences
    except Exception as e:
        logger.error(f"Error loading FASTA file: {e}")
        raise

def load_uniprot_data(csv_file: str) -> Dict[str, str]:
    """Load UniProt sequences from CSV file."""
    try:
        df = pd.read_csv(csv_file)
        accession_col = next((col for col in ['accession', 'Accession', 'ACCESSION', 'uniprot_id', 'UniProt_ID'] if col in df.columns), None)
        if accession_col is None:
            logger.error("Could not find accession column in CSV file")
            logger.info(f"Available columns: {list(df.columns)}")
            raise ValueError("Accession column not found")
        sequence_col = next((col for col in ['sequence', 'Sequence', 'SEQUENCE', 'protein_sequence', 'Protein_Sequence'] if col in df.columns), None)
        if sequence_col is None:
            logger.error("Could not find sequence column in CSV file")
            logger.info(f"Available columns: {list(df.columns)}")
            raise ValueError("Sequence column not found")
        uniprot_data = {}
        for _, row in df.iterrows():
            accession = str(row[accession_col]).strip()
            sequence = str(row[sequence_col]).strip()
            if pd.notna(accession) and pd.notna(sequence) and sequence != 'nan':
                uniprot_data[accession] = sequence
        logger.info(f"Loaded {len(uniprot_data)} UniProt sequences from {csv_file}")
        return uniprot_data
    except Exception as e:
        logger.error(f"Error loading UniProt data: {e}")
        raise

def load_accession_array(array_file: str) -> np.ndarray:
    """Load accession array from file."""
    try:
        accession_array = np.load(array_file, allow_pickle=True).astype(str)
        logger.info(f"Loaded accession array with shape {accession_array.shape}")
        return accession_array
    except Exception as e:
        logger.error(f"Error loading accession array: {e}")
        raise

def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    """
    Calculate percent identity between two protein sequences using global alignment (parasail).
    Returns percent identity (0-100).
    """
    if not seq1 or not seq2:
        return 0.0
    # Use BLOSUM62, gap open 10, gap extend 1 (standard)
    result = parasail.nw_stats_scan_16(seq1, seq2, 10, 1, parasail.blosum62)
    matches = result.matches
    length = result.length
    if length == 0:
        return 0.0
    identity = (matches / length) * 100
    return identity

def perform_pairwise_alignments(fasta_sequences: List[Tuple[str, str]], 
                               accession_array: np.ndarray,
                               uniprot_data: Dict[str, str]) -> np.ndarray:
    """Perform pairwise alignments using parasail."""
    N, k = accession_array.shape
    similarity_matrix = np.zeros((N, k))
    logger.info(f"Starting fast pairwise alignments for {N} sequences against {k} accessions each")
    for i, (seq_id, fasta_seq) in tqdm(enumerate(fasta_sequences), total=N, desc="Processing sequences"):
        for j in range(k):
            accession = str(accession_array[i, j]).strip()
            if accession in uniprot_data:
                uniprot_seq = uniprot_data[accession]
                similarity = calculate_sequence_identity(fasta_seq, uniprot_seq)
                similarity_matrix[i, j] = similarity
            else:
                logger.warning(f"Accession {accession} not found in UniProt data")
                similarity_matrix[i, j] = 0.0
    return similarity_matrix

def save_results(similarity_matrix: np.ndarray, output_file: str):
    """Save similarity matrix to file and CSV."""
    try:
        np.save(output_file, similarity_matrix)
        logger.info(f"Saved similarity matrix to {output_file}")
        csv_file = output_file.replace('.npy', '.csv')
        pd.DataFrame(similarity_matrix).to_csv(csv_file, index=False, header=False)
        logger.info(f"Saved similarity matrix as CSV to {csv_file}")
    except Exception as e:
        logger.error(f"Error saving results: {e}")
        raise

def print_summary(similarity_matrix: np.ndarray, fasta_sequences: List[Tuple[str, str]]):
    """Print summary statistics of the alignment results."""
    print("\n" + "="*60)
    print("PAIRWISE ALIGNMENT SUMMARY (FAST)")
    print("="*60)
    print(f"Number of query sequences: {len(fasta_sequences)}")
    print(f"Number of target accessions per sequence: {similarity_matrix.shape[1]}")
    print(f"Total alignments performed: {similarity_matrix.size}")
    print(f"\nSimilarity Statistics:")
    print(f"  Mean identity: {np.mean(similarity_matrix):.2f}%")
    print(f"  Median identity: {np.median(similarity_matrix):.2f}%")
    print(f"  Min identity: {np.min(similarity_matrix):.2f}%")
    print(f"  Max identity: {np.max(similarity_matrix):.2f}%")
    print(f"  Standard deviation: {np.std(similarity_matrix):.2f}%")
    print(f"\nTop matches for each sequence:")
    for i, (seq_id, _) in enumerate(fasta_sequences):
        similarities = similarity_matrix[i]
        max_idx = np.argmax(similarities)
        max_similarity = similarities[max_idx]
        print(f"  {seq_id}: {max_similarity:.2f}% (position {max_idx})")

def main():
    """Main function to run the fast pairwise alignment analysis."""
    parser = argparse.ArgumentParser(description='Perform fast pairwise sequence alignment using parasail')
    parser.add_argument('fasta_file', help='Input FASTA file with protein sequences')
    parser.add_argument('accession_array', help='NumPy array file with UniProt accessions')
    parser.add_argument('uniprot_csv', help='CSV file with UniProt sequences')
    parser.add_argument('-o', '--output', default='similarity_matrix.npy', 
                       help='Output file for similarity matrix (default: similarity_matrix.npy)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    try:
        logger.info("Loading input data...")
        fasta_sequences = load_fasta_sequences(args.fasta_file)
        accession_array = load_accession_array(args.accession_array)
        uniprot_data = load_uniprot_data(args.uniprot_csv)
        if len(fasta_sequences) != accession_array.shape[0]:
            raise ValueError(f"Number of FASTA sequences ({len(fasta_sequences)}) "
                           f"doesn't match first dimension of accession array ({accession_array.shape[0]})")
        logger.info("Performing fast pairwise alignments...")
        similarity_matrix = perform_pairwise_alignments(fasta_sequences, accession_array, uniprot_data)
        save_results(similarity_matrix, args.output)
        print_summary(similarity_matrix, fasta_sequences)
        logger.info("Analysis completed successfully!")
    except Exception as e:
        logger.error(f"Error in main execution: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 