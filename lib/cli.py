import argparse
import sys
import os
import glob
from .download_database import download_gdrive_folder
from .unnotate import unnotate

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

    args = parser.parse_args()

    if args.command == "download_db":
        download_gdrive_folder(dest_path=args.dest)
    elif args.command == "unnot":
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
            fasta_file=args.fasta_file,
            database_dir=args.database_dir,
            k=args.k,
            metric=args.metric,
            output_dir=args.output_dir,
            prefix=args.prefix
        )
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main() 
