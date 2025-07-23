import os
import zipfile
from os.path import join
import numpy as np
import pandas as pd
import logging
from .lovis4u import generate_lovis_plot

logger = logging.getLogger(__name__)


def append_product_to_gff(gff_file, full_name_csv, pval_csv=None, pval_threshold=0.01, score_csv=None, score_threshold=None):
    """
    Append product names to GFF file, optionally filtering by p-values or score.
    Args:
        gff_file: Path to GFF file
        full_name_csv: Path to full name CSV
        pval_csv: Optional path to p-value CSV
        pval_threshold: P-value threshold
        score_csv: Optional path to score CSV
        score_threshold: Score threshold (lower bound)
    """
    if not os.path.exists(full_name_csv):
        logger.warning(f"Full name CSV {full_name_csv} not found. Skipping product annotation.")
        return

    logger.info(f"Appending products to {gff_file} using {full_name_csv}")
    # Load data
    try:
        full_names = pd.read_csv(full_name_csv).iloc[:, 0].tolist()
        significant_mask = None
        if pval_csv and os.path.exists(pval_csv):
            pvals = pd.read_csv(pval_csv).values
            significant_mask = pvals < pval_threshold
            logger.info(f"Found {np.sum(significant_mask)} significant hits ({np.sum(significant_mask)/pvals.size*100:.2f}%)")
        elif score_csv and os.path.exists(score_csv) and score_threshold is not None:
            scores = pd.read_csv(score_csv).values
            significant_mask = scores >= score_threshold
            logger.info(f"Found {np.sum(significant_mask)} high-score hits ({np.sum(significant_mask)/scores.size*100:.2f}%)")
    except Exception as e:
        logger.error(f"Error reading input files: {e}")
        return

    # Process GFF file
    try:
        with open(gff_file) as f:
            lines = f.readlines()
        out_lines = []
        cds_index = 0
        for line in lines:
            if not line.strip() or line.startswith(('#', '>')):
                out_lines.append(line)
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS' and 'product=' not in parts[8]:
                # Add product if significant or no filtering
                if (significant_mask is None or 
                    (cds_index < significant_mask.shape[0] and significant_mask[cds_index, 0])):
                    product = full_names[cds_index] if cds_index < len(full_names) else ''
                    parts[8] += f';product={product}'
                else:
                    parts[8] += ';product=hypothetical protein'
                out_lines.append('\t'.join(parts) + '\n')
                cds_index += 1
            else:
                out_lines.append(line)
        with open(gff_file, 'w') as f:
            f.writelines(out_lines)
        logger.info(f"Successfully updated {gff_file}")
    except Exception as e:
        logger.error(f"Error modifying GFF: {e}")


def create_filtered_full_names_csv(full_name_csv, pval_csv, output_csv, pval_threshold=0.01):
    """
    Create a CSV file with neighbor_1 full names filtered by p-value threshold.
    
    Args:
        full_name_csv: Path to the full_name CSV file
        pval_csv: Path to the p-value CSV file
        output_csv: Path for the output filtered CSV file
        pval_threshold: P-value threshold for filtering (default: 0.01)
    """
    if not os.path.exists(full_name_csv):
        logger.warning(f"Full name CSV {full_name_csv} not found. Skipping filtered CSV creation.")
        return
    
    if not os.path.exists(pval_csv):
        logger.warning(f"P-value CSV {pval_csv} not found. Skipping filtered CSV creation.")
        return
    
    logger.info(f"Creating filtered full names CSV: {output_csv}")
    
    try:
        # Load data
        full_names_df = pd.read_csv(full_name_csv)
        pvals_df = pd.read_csv(pval_csv)
        
        # Get neighbor_1 full names and p-values
        neighbor1_full_names = full_names_df.iloc[:, 0].values  # First column
        neighbor1_pvals = pvals_df.iloc[:, 0].values  # First column
        
        # Create filtered full names (use 'hypothetical protein' for non-significant hits)
        filtered_full_names = []
        for i, (full_name, pval) in enumerate(zip(neighbor1_full_names, neighbor1_pvals)):
            filtered_full_names.append(full_name if full_name else 'hypothetical protein')
        
        # Create output DataFrame
        output_df = pd.DataFrame({
            'neighbor_1_full_name': filtered_full_names,
            'neighbor_1_pval': neighbor1_pvals,
            'is_significant': neighbor1_pvals < pval_threshold
        })
        
        # Save to CSV
        output_df.to_csv(output_csv, index=False)
        
        significant_count = np.sum(neighbor1_pvals < pval_threshold)
        total_count = len(neighbor1_pvals)
        logger.info(f"Created filtered CSV with {significant_count}/{total_count} significant hits ({significant_count/total_count*100:.2f}%)")
        
    except Exception as e:
        logger.error(f"Error creating filtered CSV: {e}")


def generate_annotation_pdf(gff_file, output_pdf_path):
    """
    Generate lovis4u PDF annotation visualization.
    
    Args:
        gff_file: Path to the GFF file
        output_pdf_path: Path for the output PDF file
    """
    if not os.path.exists(gff_file):
        logger.warning(f"GFF file {gff_file} not found. Skipping PDF generation.")
        return
    
    logger.info(f"Generating lovis4u PDF at {output_pdf_path} from {gff_file}")
    
    try:
        generate_lovis_plot(gff_path=gff_file, output_pdf_path=output_pdf_path)
        logger.info(f"Successfully generated annotation PDF: {output_pdf_path}")
    except Exception as e:
        logger.error(f"Error generating PDF: {e}")


def summarize_outputs(output_dir, prefix, input_type="protein", gff_prefix=None, pval_csv=None, pval_threshold=0.01, score_csv=None, score_threshold=None):
    """
    Perform post-processing tasks for Unnotate outputs.
    Args:
        output_dir: Directory containing the output files
        prefix: Prefix used for output files
        input_type: Type of input ("dna" or "protein")
        gff_prefix: Prefix for GFF file (if different from main prefix, e.g., for DNA input)
        pval_csv: Optional path to p-value CSV file for filtering products
        pval_threshold: P-value threshold for including products (default: 0.01)
        score_csv: Optional path to score CSV file for filtering products
        score_threshold: Score threshold (lower bound) for including products
    """
    logger.info(f"Starting post-processing for {input_type} input with prefix '{prefix}'")
    # Define file paths
    if gff_prefix:
        gff_file = os.path.join(output_dir, f"{gff_prefix}.gff")
    else:
        gff_file = os.path.join(output_dir, f"{prefix}.gff")
    full_name_csv = os.path.join(output_dir, f"{prefix}_full_name.csv")
    output_pdf_path = os.path.join(output_dir, f"{prefix}_annotation.pdf")
    filtered_csv_path = os.path.join(output_dir, f"{prefix}_filtered_full_names.csv")
    # Filtering logic
    filter_type = None
    filter_matrix = None
    if pval_csv and os.path.exists(pval_csv):
        logger.info(f"Using p-value CSV {pval_csv} with threshold {pval_threshold}")
        filter_type = 'pval'
        filter_matrix = pd.read_csv(pval_csv).values < pval_threshold
    elif score_csv and os.path.exists(score_csv) and score_threshold is not None:
        logger.info(f"Using score CSV {score_csv} with threshold {score_threshold}")
        filter_type = 'score'
        filter_matrix = pd.read_csv(score_csv).values >= score_threshold
    else:
        logger.info("No p-value or score filtering applied.")
    # Create filtered full names CSV if filtering is applied
    if filter_type and os.path.exists(full_name_csv):
        full_names = pd.read_csv(full_name_csv).values
        filtered_full_names = np.where(filter_matrix, full_names, 'hypothetical protein')
        # Save filtered CSV
        pd.DataFrame(filtered_full_names, columns=[f"neighbor_{i+1}" for i in range(filtered_full_names.shape[1])]).to_csv(filtered_csv_path, index=False)
        logger.info(f"Created filtered CSV at {filtered_csv_path}")
    # For DNA input, handle GFF and PDF generation
    if input_type == "dna":
        # Append product information to GFF
        if os.path.exists(gff_file):
            append_product_to_gff(gff_file, full_name_csv, pval_csv if filter_type == 'pval' else None, pval_threshold, score_csv, score_threshold)
        else:
            logger.warning(f"GFF file {gff_file} not found for DNA input")
        # Generate lovis4u PDF
        if os.path.exists(gff_file):
            generate_annotation_pdf(gff_file, output_pdf_path)
        else:
            logger.warning(f"Cannot generate PDF: GFF file {gff_file} not found")
    else:
        logger.info("Protein input detected. No additional post-processing required.")
    logger.info("Post-processing completed")


def save_results(output_dir, df, accessions, similarity, percent_identity_matrix, score_matrix, prefix="unnotated"):
    """
    Save results to output directory and create zip file. All arrays are saved in a single .npz file.
    Also outputs full_name and domain matrices corresponding to the accession matrix.
    """
    os.makedirs(output_dir, exist_ok=True)
    # Save individual files
    csv_path = join(output_dir, f"{prefix}_uniprot.csv")
    npz_path = join(output_dir, f"{prefix}_results.npz")
    # Create CSV files for sequence identity and cosine similarity
    similarity_csv_path = join(output_dir, f"{prefix}_cosine_similarity.csv")
    identity_csv_path = join(output_dir, f"{prefix}_sequence_identity.csv")
    score_csv_path = join(output_dir, f"{prefix}_score.csv")
    accession_csv_path = join(output_dir, f"{prefix}_accession.csv")
    full_name_csv_path = join(output_dir, f"{prefix}_full_name.csv")
    domain_csv_path = join(output_dir, f"{prefix}_domain.csv")
    sequence_length_csv_path = join(output_dir, f"{prefix}_sequence_length.csv")
    df.to_csv(csv_path, index=False)
    # Save all arrays to npz, including score if provided
    np.savez(npz_path, accession=accessions, cosine_similarity=similarity, sequence_identity=percent_identity_matrix, score=score_matrix)
    # Save cosine similarity as CSV
    similarity_df = pd.DataFrame(similarity, columns=[f"neighbor_{i+1}" for i in range(similarity.shape[1])])
    similarity_df.to_csv(similarity_csv_path, index=False)
    # Save sequence identity as CSV
    identity_df = pd.DataFrame(percent_identity_matrix, columns=[f"neighbor_{i+1}" for i in range(percent_identity_matrix.shape[1])])
    identity_df.to_csv(identity_csv_path, index=False)
    # Save accession as CSV (ensure type is string)
    accession_df = pd.DataFrame(accessions.astype(str), columns=[f"neighbor_{i+1}" for i in range(accessions.shape[1])])
    accession_df.to_csv(accession_csv_path, index=False)
    # Build mapping from accession to full_name, taxonomy domain, and sequence length
    acc_to_full_name = {row["accession"]: row["full_name"] if pd.notna(row["full_name"]) else "" for _, row in df.iterrows()}
    acc_to_taxonomy = {row["accession"]: row["taxonomy_lineage"] if pd.notna(row["taxonomy_lineage"]) else "" for _, row in df.iterrows()}
    acc_to_sequence_length = {row["accession"]: (row["sequence_length"]) if pd.notna(row["sequence_length"]) else "" for _, row in df.iterrows()}
    # Build full_name and domain matrices
    full_name_matrix = np.empty_like(accessions, dtype=object)
    taxa_domain_matrix = np.empty_like(accessions, dtype=object)
    sequence_length_matrix = np.empty_like(accessions, dtype=int)
    for i in range(accessions.shape[0]):
        for j in range(accessions.shape[1]):
            acc = accessions[i, j]
            full_name_matrix[i, j] = acc_to_full_name.get(acc, "Unknown")
            sequence_length_matrix[i, j] = acc_to_sequence_length.get(acc, "")
            taxonomy = acc_to_taxonomy.get(acc, "")
            if taxonomy:
                domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                taxa_domain_matrix[i, j] = domain if domain else "Unknown"
            else:
                taxa_domain_matrix[i, j] = "Unknown"
    # Save as CSV
    full_name_df = pd.DataFrame(full_name_matrix, columns=[f"neighbor_{i+1}" for i in range(full_name_matrix.shape[1])])
    full_name_df.to_csv(full_name_csv_path, index=False)
    domain_df = pd.DataFrame(taxa_domain_matrix, columns=[f"neighbor_{i+1}" for i in range(taxa_domain_matrix.shape[1])])
    domain_df.to_csv(domain_csv_path, index=False)
    sequence_length_df = pd.DataFrame(sequence_length_matrix, columns=[f"neighbor_{i+1}" for i in range(sequence_length_matrix.shape[1])])
    sequence_length_df.to_csv(sequence_length_csv_path, index=False)
    score_df = pd.DataFrame(score_matrix, columns=[f"neighbor_{i+1}" for i in range(score_matrix.shape[1])])
    score_df.to_csv(score_csv_path, index=False)
    # Create zip file containing only uniprot.csv and results.npz (not the individual CSV files)
    zip_path = join(output_dir, f"{prefix}_streamlit.zip")
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        zipf.write(csv_path, arcname=f"{prefix}_uniprot.csv")
        zipf.write(npz_path, arcname=f"{prefix}_results.npz")
    logger.info(f"Saved UniProt data in {csv_path}")
    logger.info(f"Saved cosine similarity matrix in {similarity_csv_path}")
    logger.info(f"Saved sequence identity matrix in {identity_csv_path}")
    logger.info(f"Saved accession matrix in {accession_csv_path}")
    logger.info(f"Saved full_name matrix in {full_name_csv_path}")
    logger.info(f"Saved domain matrix in {domain_csv_path}")
    logger.info(f"Saved sequence length matrix in {sequence_length_csv_path}")
    logger.info(f"Saved all results in {npz_path}")
    logger.info(f"Created zip archive of Streamlit-compatible files at {zip_path}")


def main():
    """Command line interface for summarize outputs."""
    import argparse
    parser = argparse.ArgumentParser(description="Post-process Unnotate outputs")
    parser.add_argument("--output-dir", required=True, help="Directory containing output files")
    parser.add_argument("--prefix", required=True, help="Prefix for output files")
    parser.add_argument("--input-type", choices=["protein", "dna"], default="protein", help="Type of input (protein or dna)")
    parser.add_argument("--gff-prefix", default=None, help="Prefix for GFF file (optional)")
    parser.add_argument("--pval-csv", default=None, help="Path to p-value CSV file (optional)")
    parser.add_argument("--pval-threshold", type=float, default=0.01, help="P-value threshold (default: 0.01)")
    parser.add_argument("--score-csv", default=None, help="Path to score CSV file (optional)")
    parser.add_argument("--score-threshold", type=float, default=None, help="Score threshold (lower bound, optional)")
    args = parser.parse_args()
    try:
        summarize_outputs(
            output_dir=args.output_dir,
            prefix=args.prefix,
            input_type=args.input_type,
            gff_prefix=args.gff_prefix,
            pval_csv=args.pval_csv,
            pval_threshold=args.pval_threshold,
            score_csv=args.score_csv,
            score_threshold=args.score_threshold
        )
    except Exception as e:
        logger.error(f"Error in post-processing: {e}")
        raise


if __name__ == '__main__':
    main() 