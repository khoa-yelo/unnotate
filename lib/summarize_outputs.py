import os
import numpy as np
import pandas as pd
import logging
from .lovis4u import generate_lovis_plot

logger = logging.getLogger(__name__)


def append_product_to_gff(gff_file, full_name_csv, pval_csv=None, pval_threshold=0.01):
    """Append product names to GFF file, optionally filtering by p-values."""
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
                # Add product if significant or no p-value filtering
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


def summarize_outputs(output_dir, prefix, input_type="protein", gff_prefix=None, pval_csv=None, pval_threshold=0.01):
    """
    Perform post-processing tasks for Unnotate outputs.
    
    Args:
        output_dir: Directory containing the output files
        prefix: Prefix used for output files
        input_type: Type of input ("dna" or "protein")
        gff_prefix: Prefix for GFF file (if different from main prefix, e.g., for DNA input)
        pval_csv: Optional path to p-value CSV file for filtering products
        pval_threshold: P-value threshold for including products (default: 0.01)
    """
    logger.info(f"Starting post-processing for {input_type} input with prefix '{prefix}'")
    
    # Define file paths
    if gff_prefix:
        gff_file = os.path.join(output_dir, f"{gff_prefix}.gff")
    else:
        gff_file = os.path.join(output_dir, f"{prefix}.gff")
    
    full_name_csv = os.path.join(output_dir, f"{prefix}_full_name.csv")
    output_pdf_path = os.path.join(output_dir, f"{prefix}_annotation.pdf")
    
    # For DNA input, we need to handle GFF and PDF generation
    if input_type == "dna":
        # Append product information to GFF
        if os.path.exists(gff_file):
            append_product_to_gff(gff_file, full_name_csv, pval_csv, pval_threshold)
        else:
            logger.warning(f"GFF file {gff_file} not found for DNA input")
        
        # Generate lovis4u PDF
        if os.path.exists(gff_file):
            generate_annotation_pdf(gff_file, output_pdf_path)
        else:
            logger.warning(f"Cannot generate PDF: GFF file {gff_file} not found")
    
    # For protein input, no additional processing is needed
    else:
        logger.info("Protein input detected. No additional post-processing required.")
    
    logger.info("Post-processing completed")


def main():
    """Command line interface for summarize outputs."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Post-process Unnotate outputs")
    parser.add_argument("--output-dir", required=True, help="Directory containing output files")
    parser.add_argument("--prefix", required=True, help="Prefix for output files")
    parser.add_argument("--input-type", choices=["protein", "dna"], default="protein", 
                       help="Type of input (protein or dna)")
    
    args = parser.parse_args()
    
    try:
        summarize_outputs(
            output_dir=args.output_dir,
            prefix=args.prefix,
            input_type=args.input_type
        )
    except Exception as e:
        logger.error(f"Error in post-processing: {e}")
        raise


if __name__ == '__main__':
    main() 