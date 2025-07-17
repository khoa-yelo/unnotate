import numpy as np
import pandas as pd
import os
import logging
from .pval import upper_p_value

logger = logging.getLogger(__name__)

def load_random_background(database_dir):
    """Load random background similarities from random_proteins.npz."""
    npz_path = os.path.join(database_dir, 'random_proteins.npz')
    if not os.path.exists(npz_path):
        raise FileNotFoundError(f"random_proteins.npz not found in {database_dir}")
    random_out = np.load(npz_path)
    x, y = random_out["query_lengths"], random_out["distances"][:,0]
    return x, y


def compute_summary_statistics(database_dir, cosine_similarity_csv, sequence_length_csv=None, output_dir=None, prefix=None):
    """
    Compute summary statistics and p-values for Unnotate results.
    
    Args:
        database_dir: Directory containing random_proteins.npz
        cosine_similarity_csv: Path to cosine similarity CSV file
        sequence_length_csv: Optional path to sequence length CSV file
        output_file: Optional output file path
    
    Returns:
        DataFrame with summary statistics
    """
    query_sim = pd.read_csv(cosine_similarity_csv).values
    query_length = pd.read_csv(sequence_length_csv).values
    x, y = load_random_background(database_dir)

    logger.info("Computing p-values...")
    pvals = []
    for i in range(len(query_length)):
        for j in range(len(query_sim[i])):
            pvals.append(upper_p_value(query_length[i][j], query_sim[i][j], x, y))
    pvals = np.array(pvals).reshape(len(query_length), len(query_sim[0]))
    
    # Create pval DataFrame with pvals
    df_out = pd.DataFrame(pvals)
    # Print summary
    logger.info(f"Summary statistics computed for {len(query_length)} queries")
    
    # Save to file if requested
    if output_dir:
        df_out.columns = [f"neighbor_{i+1}" for i in range(len(df_out.columns))]
        df_out.to_csv(os.path.join(output_dir, f"{prefix}_pvals.csv"), index=False)
        logger.info(f"Summary statistics written to {os.path.join(output_dir, f'{prefix}_pvals.csv')}")
    
    return df_out


def main():
    """Command line interface for summary statistics."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Compute summary statistics and p-values for Unnotate results.")
    parser.add_argument('--database-dir', type=str, required=True, help='Directory containing random_proteins.npz')
    parser.add_argument('--cosine-similarity-csv', type=str, required=True, help='CSV file with cosine similarity matrix (N×k)')
    parser.add_argument('--sequence-length-csv', type=str, default=None, help='CSV file with sequence lengths (optional)')
    parser.add_argument('--output-dir', type=str, default=None, help='Output directory for summary statistics (optional)')
    
    args = parser.parse_args()
    
    try:
        compute_summary_statistics(
            database_dir=args.database_dir,
            cosine_similarity_csv=args.cosine_similarity_csv,
            sequence_length_csv=args.sequence_length_csv,
            output_dir=args.output_dir,
            prefix=args.prefix
        )
        
    except Exception as e:
        logger.error(f"Error computing summary statistics: {e}")
        raise


if __name__ == '__main__':
    main()
