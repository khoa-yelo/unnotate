#!/usr/bin/env python3
"""
Pairwise Sequence Alignment Script

This script performs pairwise sequence alignment between protein sequences from a FASTA file
and sequences from UniProt accessions, reporting percent similarity.

Input:
- FASTA file containing N protein sequences
- NumPy array accession_array with shape (N, k) containing UniProt accessions
- CSV file with UniProt sequences (parsed_uniprot_swiss_data.csv)

Output:
- Percent similarity matrix for each protein sequence against its k accessions
"""

import numpy as np
import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
import argparse
import sys
from typing import List, Tuple, Dict
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_fasta_sequences(fasta_file: str) -> List[Tuple[str, str]]:
    """
    Load protein sequences from FASTA file.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        List of tuples (sequence_id, sequence)
    """
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
    """
    Load UniProt sequences from CSV file.
    
    Args:
        csv_file: Path to CSV file containing UniProt data
        
    Returns:
        Dictionary mapping accession to sequence
    """
    try:
        # Try different possible column names for accession and sequence
        df = pd.read_csv(csv_file)
        
        # Look for accession column (common names)
        accession_col = None
        for col in ['accession', 'Accession', 'ACCESSION', 'uniprot_id', 'UniProt_ID']:
            if col in df.columns:
                accession_col = col
                break
        
        if accession_col is None:
            logger.error("Could not find accession column in CSV file")
            logger.info(f"Available columns: {list(df.columns)}")
            raise ValueError("Accession column not found")
        
        # Look for sequence column (common names)
        sequence_col = None
        for col in ['sequence', 'Sequence', 'SEQUENCE', 'protein_sequence', 'Protein_Sequence']:
            if col in df.columns:
                sequence_col = col
                break
        
        if sequence_col is None:
            logger.error("Could not find sequence column in CSV file")
            logger.info(f"Available columns: {list(df.columns)}")
            raise ValueError("Sequence column not found")
        
        # Create dictionary
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
    """
    Load accession array from file.
    
    Args:
        array_file: Path to NumPy array file
        
    Returns:
        NumPy array with shape (N, k)
    """
    try:
        accession_array = np.load(array_file, allow_pickle=True).astype(str)
        logger.info(f"Loaded accession array with shape {accession_array.shape}")
        return accession_array
    except Exception as e:
        logger.error(f"Error loading accession array: {e}")
        raise


def calculate_sequence_similarity(seq1: str, seq2: str) -> float:
    """
    Calculate percent similarity between two protein sequences using global alignment.
    
    Args:
        seq1: First protein sequence
        seq2: Second protein sequence
        
    Returns:
        Percent similarity (0-100)
    """
    if not seq1 or not seq2:
        return 0.0
    
    # Perform global alignment
    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
    
    if not alignments:
        return 0.0
    
    # Get the best alignment
    best_alignment = alignments[0]
    aligned_seq1, aligned_seq2 = best_alignment[0], best_alignment[1]
    
    # Calculate similarity
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
    total_length = len(aligned_seq1)
    
    if total_length == 0:
        return 0.0
    
    similarity = (matches / total_length) * 100
    return similarity


def perform_pairwise_alignments(fasta_sequences: List[Tuple[str, str]], 
                               accession_array: np.ndarray,
                               uniprot_data: Dict[str, str]) -> np.ndarray:
    """
    Perform pairwise alignments between FASTA sequences and UniProt sequences.
    
    Args:
        fasta_sequences: List of (sequence_id, sequence) tuples
        accession_array: NumPy array with shape (N, k) containing accessions
        uniprot_data: Dictionary mapping accession to sequence
        
    Returns:
        NumPy array with similarity scores, shape (N, k)
    """
    N, k = accession_array.shape
    similarity_matrix = np.zeros((N, k))
    
    logger.info(f"Starting pairwise alignments for {N} sequences against {k} accessions each")
    
    for i, (seq_id, fasta_seq) in enumerate(fasta_sequences):
        logger.info(f"Processing sequence {i+1}/{N}: {seq_id}")
        
        for j in range(k):
            accession = str(accession_array[i, j]).strip()
            
            if accession in uniprot_data:
                uniprot_seq = uniprot_data[accession]
                similarity = calculate_sequence_similarity(fasta_seq, uniprot_seq)
                similarity_matrix[i, j] = similarity
            else:
                logger.warning(f"Accession {accession} not found in UniProt data")
                similarity_matrix[i, j] = 0.0
    
    return similarity_matrix


def save_results(similarity_matrix: np.ndarray, output_file: str):
    """
    Save similarity matrix to file.
    
    Args:
        similarity_matrix: NumPy array with similarity scores
        output_file: Path to output file
    """
    try:
        np.save(output_file, similarity_matrix)
        logger.info(f"Saved similarity matrix to {output_file}")
        
        # Also save as CSV for easier viewing
        csv_file = output_file.replace('.npy', '.csv')
        df = pd.DataFrame(similarity_matrix)
        df.to_csv(csv_file, index=False, header=False)
        logger.info(f"Saved similarity matrix as CSV to {csv_file}")
        
    except Exception as e:
        logger.error(f"Error saving results: {e}")
        raise


def print_summary(similarity_matrix: np.ndarray, fasta_sequences: List[Tuple[str, str]]):
    """
    Print summary statistics of the alignment results.
    
    Args:
        similarity_matrix: NumPy array with similarity scores
        fasta_sequences: List of (sequence_id, sequence) tuples
    """
    print("\n" + "="*60)
    print("PAIRWISE ALIGNMENT SUMMARY")
    print("="*60)
    
    print(f"Number of query sequences: {len(fasta_sequences)}")
    print(f"Number of target accessions per sequence: {similarity_matrix.shape[1]}")
    print(f"Total alignments performed: {similarity_matrix.size}")
    
    print(f"\nSimilarity Statistics:")
    print(f"  Mean similarity: {np.mean(similarity_matrix):.2f}%")
    print(f"  Median similarity: {np.median(similarity_matrix):.2f}%")
    print(f"  Min similarity: {np.min(similarity_matrix):.2f}%")
    print(f"  Max similarity: {np.max(similarity_matrix):.2f}%")
    print(f"  Standard deviation: {np.std(similarity_matrix):.2f}%")
    
    # Show top matches for each sequence
    print(f"\nTop matches for each sequence:")
    for i, (seq_id, _) in enumerate(fasta_sequences):
        similarities = similarity_matrix[i]
        max_idx = np.argmax(similarities)
        max_similarity = similarities[max_idx]
        print(f"  {seq_id}: {max_similarity:.2f}% (position {max_idx})")


def main():
    """Main function to run the pairwise alignment analysis."""
    parser = argparse.ArgumentParser(description='Perform pairwise sequence alignment')
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
        # Load data
        logger.info("Loading input data...")
        fasta_sequences = load_fasta_sequences(args.fasta_file)
        accession_array = load_accession_array(args.accession_array)
        uniprot_data = load_uniprot_data(args.uniprot_csv)
        
        # Validate input dimensions
        if len(fasta_sequences) != accession_array.shape[0]:
            raise ValueError(f"Number of FASTA sequences ({len(fasta_sequences)}) "
                           f"doesn't match first dimension of accession array ({accession_array.shape[0]})")
        
        # Perform alignments
        logger.info("Performing pairwise alignments...")
        similarity_matrix = perform_pairwise_alignments(fasta_sequences, accession_array, uniprot_data)
        
        # Save results
        save_results(similarity_matrix, args.output)
        
        # Print summary
        print_summary(similarity_matrix, fasta_sequences)
        
        logger.info("Analysis completed successfully!")
        
    except Exception as e:
        logger.error(f"Error in main execution: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 