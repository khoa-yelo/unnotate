import argparse
import sys
import os
import glob
import logging
import torch
import pathlib
import tempfile
from Bio import SeqIO
from .download_database import download_gdrive_folder
from .unnotate import unnotate

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
            # Find the .faa and .gff files
            prefix = os.path.splitext(os.path.basename(fasta_file))[0]
            faa_file = os.path.join(args.output_dir, f"{prefix}_proteins.faa")
            gff_file = os.path.join(args.output_dir, f"{prefix}.gff")
            if not os.path.exists(faa_file):
                logger.error(f"Protein FASTA file {faa_file} not found after translation.")
                sys.exit(1)
            fasta_file = faa_file
            logger.info(f"Using translated protein FASTA: {fasta_file}")
            # Only for DNA input: generate lovis4u PDF
            lovis4u_pdf = os.path.join(args.output_dir, f"{prefix}_annotation.pdf")
            logger.info(f"Generating lovis4u PDF at {lovis4u_pdf} from {gff_file}")
            generate_lovis_plot(gff_path=gff_file, output_pdf_path=lovis4u_pdf)
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
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main() 
