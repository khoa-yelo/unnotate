import argparse
import sys
from .download_database import download_gdrive_folder
from .unnotate import unnotate

def main():
    parser = argparse.ArgumentParser(prog="unnotate")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # download_db subcommand
    parser_db = subparsers.add_parser("download_db", help="Download the database from Google Drive")
    parser_db.add_argument("--dest", default="./database", help="Destination directory")

    # annot subcommand
    parser_annot = subparsers.add_parser("annot", help="Annotate protein sequences")
    parser_annot.add_argument("--fasta_file", required=True)
    parser_annot.add_argument("--embeddings_db", required=True)
    parser_annot.add_argument("--uniprot_db", required=True)
    parser_annot.add_argument("--k", type=int, default=20)
    parser_annot.add_argument("--metric", default="mean_middle_layer_12")
    parser_annot.add_argument("--output_dir", required=True)
    parser_annot.add_argument("--prefix", default="unnotated")

    args = parser.parse_args()

    if args.command == "download_db":
        download_gdrive_folder(dest_path=args.dest)
    elif args.command == "unnot":
        unnotate(
            fasta_file=args.fasta_file,
            embeddings_db=args.embeddings_db,
            uniprot_db=args.uniprot_db,
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
