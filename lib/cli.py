import argparse
import sys
import os
import glob
import logging
import torch
import pandas as pd
import pathlib
from Bio import SeqIO
from .download_database import download_gdrive_folder
from .unnotate import unnotate
from .summary_statistics import compute_summary_statistics

# Local import for DNA translation (pyrodigal)
from .pyrodigal import run_pyrodigal_gv
from .lovis4u import generate_lovis_plot

def detect_input_type(fasta_file, user_input_type):
    """
    Detect if the input FASTA is DNA or protein.
    If user_input_type is not 'auto', return it.
    Otherwise, read the first sequence and check vocabulary.
    """
    if user_input_type != "auto":
        return user_input_type
    with open(fasta_file) as handle:
        first_record = next(SeqIO.parse(handle, "fasta"))
        seq = str(first_record.seq).upper()
        dna_letters = set("ATCGN")
        seq_letters = set(seq)
        # If all letters are in DNA set, treat as DNA, else protein
        if seq_letters <= dna_letters:
            return "dna"
        else:
            return "protein"

def append_product_to_gff(gff_file, full_name_csv):
    """Append product=[full_name] to each CDS in GFF using first column of full_name CSV."""
    full_names = pd.read_csv(full_name_csv).iloc[:, 0].tolist()
    with open(gff_file) as f:
        lines = f.readlines()
    out, idx = [], 0
    for line in lines:
        if not (line.strip() and not line.startswith(('#', '>'))):
            out.append(line)
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 9 and parts[2] == 'CDS':
            if 'product=' not in parts[8]:
                product = full_names[idx] if idx < len(full_names) else ''
                parts[8] += f';product={product}'
            out.append('\t'.join(parts) + '\n')
            idx += 1
        else:
            out.append(line)
    with open(gff_file, 'w') as f:
        f.writelines(out)

def main():
    parser = argparse.ArgumentParser(prog="unnotate")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # download_db subcommand
    parser_db = subparsers.add_parser("download_db", help="Download the database from Google Drive")
    parser_db.add_argument("--dest", default="./database", help="Destination directory")

    # unnot subcommand
    parser_annot = subparsers.add_parser("unnot", help="Annotate protein sequences")
    parser_annot.add_argument("--fasta-file", required=True)
    parser_annot.add_argument("--database-dir", required=True, help="Directory containing the embeddings HDF5 and UniProt CSV files")
    parser_annot.add_argument("--k", type=int, default=20)
    parser_annot.add_argument("--metric", default="mean_middle_layer_12")
    parser_annot.add_argument("--output-dir", required=True)
    parser_annot.add_argument("--prefix", default="unnotated")
    parser_annot.add_argument("--faiss-metric", default="cosine")
    parser_annot.add_argument("--cpu", action="store_true")
    parser_annot.add_argument("--input-type", choices=["auto", "protein", "dna"], default="auto", help="Type of input: protein (amino acid FASTA), dna (nucleotide FASTA), or auto (infer from sequence)")

    args = parser.parse_args()

    logger = logging.getLogger(__name__)
    logger.info(f"Starting Unnotate with args: {args}")
    logger.info(f"CUDA available: {torch.cuda.is_available()}")
    logger.info(f"Using GPU: {not args.cpu}")
    if args.cpu:
        use_gpu = False
    else:
        use_gpu = torch.cuda.is_available()
        if not use_gpu:
            logger.warning("GPU is not available, using CPU instead")
    if torch.cuda.is_available() and not use_gpu:
        logger.warning("CUDA is available, but use_gpu is not set. Try running without --cpu to use GPU - much faster.")

    if args.command == "download_db":
        download_gdrive_folder(dest_path=args.dest)
    elif args.command == "unnot":
        fasta_file = args.fasta_file
        # Detect input type (DNA or protein)
        input_type = detect_input_type(fasta_file, args.input_type)
        logger.info(f"Detected input type: {input_type}")
        if input_type == "dna":
            # Translate DNA to protein using pyrodigal, output to args.output_dir
            logger.info(f"Translating DNA FASTA {fasta_file} to protein using pyrodigal_gv...")
            run_pyrodigal_gv(fasta_path=fasta_file, output_dir=pathlib.Path(args.output_dir))
            fasta_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
            # Use consistent prefix for all outputs
            faa_file = os.path.join(args.output_dir, f"{fasta_prefix}_proteins.faa")
            gff_file = os.path.join(args.output_dir, f"{fasta_prefix}.gff")
            full_name_csv = os.path.join(args.output_dir, f"{args.prefix}_full_name.csv")
            lovis4u_pdf = os.path.join(args.output_dir, f"{args.prefix}_annotation.pdf")
            # Check for expected outputs
            if not os.path.exists(faa_file):
                logger.error(f"Protein FASTA file {faa_file} not found after translation.")
                sys.exit(1)
            fasta_file = faa_file
            logger.info(f"Using translated protein FASTA: {fasta_file}")
        # For protein input, do nothing extra (no dummy GFF, no lovis4u)
        # Find embeddings and CSV in the database dir
        h5_files = glob.glob(os.path.join(args.database_dir, "*.h5"))
        csv_files = glob.glob(os.path.join(args.database_dir, "*.csv"))
        if not h5_files:
            print(f"No .h5 embeddings file found in {args.database_dir}", file=sys.stderr)
            sys.exit(1)
        if not csv_files:
            print(f"No .csv UniProt file found in {args.database_dir}", file=sys.stderr)
            sys.exit(1)
        # Run annotation pipeline
        unnotate(
            fasta_file=fasta_file,
            database_dir=args.database_dir,
            k=args.k,
            metric=args.metric,
            output_dir=args.output_dir,
            prefix=args.prefix,
            faiss_metric=args.faiss_metric,
            use_gpu=use_gpu
        )
        if input_type == "dna":
            # Append product to GFF using full_name.csv
            if os.path.exists(full_name_csv):
                logger.info(f"Appending product=[full_name] to {gff_file} using {full_name_csv}")
                append_product_to_gff(gff_file, full_name_csv)
            # Generate lovis4u PDF
            logger.info(f"Generating lovis4u PDF at {lovis4u_pdf} from {gff_file}")
            generate_lovis_plot(gff_path=gff_file, output_pdf_path=lovis4u_pdf)
        logger.info("Computing summary statistics...")
        compute_summary_statistics(
            database_dir=args.database_dir,
            cosine_similarity_csv=os.path.join(args.output_dir, f"{args.prefix}_cosine_similarity.csv"),
            sequence_length_csv=os.path.join(args.output_dir, f"{args.prefix}_sequence_length.csv"),
            output_dir=args.output_dir,
            prefix=args.prefix
        )
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main() 
