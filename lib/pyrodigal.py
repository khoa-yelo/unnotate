import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from typing import Tuple

from Bio import SeqIO
import pyrodigal_gv


def run_pyrodigal_gv(
    fasta_path: Path,
    output_dir: Path,
    threads: int = 1,
    meta: bool = True,
    logger: logging.Logger = logging.getLogger(__name__),
) -> None:
    """
    Predict viral CDS from an input FASTA using pyrodigal_gv
    and write out GFF, nucleotide, and amino-acid FASTA files.

    Args:
        fasta_path: Path to the input FASTA file.
        output_dir: Directory where outputs will be written.
        threads: Number of worker threads.
        meta: Whether to run in metagenomic mode.
        coding_table: Genetic code table number for Prodigal.
        logger: Optional logger for debug/info messages.
    """
    fasta_path = Path(fasta_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = fasta_path.stem
    gff_file = output_dir / f"{prefix}.gff"
    nuc_fasta = output_dir / f"{prefix}_genes.fna"
    aa_fasta  = output_dir / f"{prefix}_proteins.faa"

    logger.info(f"Initializing ViralGeneFinder (meta={meta})")
    orf_finder = pyrodigal_gv.ViralGeneFinder(meta=meta)

    def _find_genes(record):
        genes = orf_finder.find_genes(str(record.seq))
        return record.id, genes

    records = SeqIO.parse(fasta_path, "fasta")

    logger.info("Starting gene finding on %s with %d threads", fasta_path, threads)
    with ThreadPoolExecutor(max_workers=threads) as pool, \
         open(gff_file, "w") as gff_out, \
         open(nuc_fasta, "w") as nuc_out, \
         open(aa_fasta, "w") as aa_out:

        for record_id, genes in pool.map(_find_genes, records):
            logger.debug("Writing results for %s: %d genes", record_id, len(genes))
            genes.write_gff(
                gff_out,
                sequence_id=record_id,
                include_translation_table=True
            )
            # append fasta of the sequence at the end of gff
            gff_out.write("##FASTA\n")
            gff_out.write(f">{record_id}\n")
            gff_out.write(f"{str(genes.sequence)}\n")
            genes.write_genes(nuc_out, sequence_id=record_id)
            genes.write_translations(aa_out, sequence_id=record_id)

    logger.info(f"Gene Finding Complete! Outputs written to {output_dir}")