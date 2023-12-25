"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from Bio.SeqIO.FastaIO import SimpleFastaParser  # type: ignore
from os.path import isfile
from tqdm import tqdm  # type: ignore
from typing import List, Tuple


def extract_kmers(input_path: str, kmer_size: int) -> List[Tuple[str, str]]:
    if not isfile(input_path):
        raise AssertionError
    with open(input_path) as FH:
        oligos_list: List[Tuple[str, str]] = []
        for record_header, record_sequence in tqdm(
            SimpleFastaParser(FH), desc="Parsing region", leave=False
        ):
            for i in tqdm(
                range(len(record_sequence) - kmer_size + 1),
                desc="Extracting oligos",
                leave=False,
            ):
                oligo_sequence = str(record_sequence)[slice(i, i + kmer_size)]
                if "N" in oligo_sequence:
                    continue
                oligos_list.append(
                    (
                        f"{record_header}|{i+1}:{i+kmer_size+1}",
                        oligo_sequence,
                    )
                )
    return oligos_list
