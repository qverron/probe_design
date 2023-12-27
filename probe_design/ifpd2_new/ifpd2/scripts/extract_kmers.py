"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from Bio.Seq import reverse_complement  # type: ignore
import click  # type: ignore
from ..const import CONTEXT_SETTINGS
from ..fasta import extract_kmers
from ..io import write_oligos
from ..oligo import select_by_GC
import logging
from os.path import isdir, isfile, join as path_join
from os.path import basename, normpath, splitext


#@click.command(
#    name="extract_kmers",
#    context_settings=CONTEXT_SETTINGS,
#    help="Generate oligonucleotides K-mers from FASTA",
#)
#@click.argument("input_path", metavar="INPUT_FASTA", type=click.Path(exists=True))
#@click.argument("output_path", metavar="OUTPUT_DIRECTORY", type=click.Path(exists=True))
#@click.argument("kmer_size", metavar="KMER_LENGTH", type=click.INT)
#@click.option('-gcfilter')
def main(input_path: str, output_path: str, kmer_size: int, gcfilter: bool) -> None:
    if not isfile(input_path):
        raise AssertionError
    if not isdir(output_path):
        raise AssertionError

    logging.info(f"Input  : {input_path}")
    logging.info(f"Output : {output_path}")
    logging.info(f"Length : {kmer_size}")

    oligos_list = extract_kmers(input_path, kmer_size)
    logging.info(f"Extracted {len(oligos_list)} sequences")

#    if (gcfilter=="False"):
    if (not gcfilter):
        logging.info(f"Skipped filtering on GC-content.")
        valid_oligos = oligos_list
    else:
        valid_oligos = select_by_GC(oligos_list, kmer_size)
        logging.info(f"{len(valid_oligos)} sequences with correct GC content")
            
    base, _ = splitext(basename(input_path))
    write_oligos(
        path_join(normpath(output_path), f"{base}.GC35to85_Reference.fa"),
        iter(valid_oligos),
    )

    write_oligos(
        path_join(normpath(output_path), f"{base}.GC35to85_RevCompl.fa"),
        ((h, reverse_complement(s)) for h, s in valid_oligos),
        "writing rc",
    )

    logging.info("Done. :thumbs_up: :smiley:")
    logging.shutdown()
