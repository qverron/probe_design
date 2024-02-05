#!/usr/bin/python3

from Bio.SeqIO.FastaIO import SimpleFastaParser
from ifpd2q.scripts.extract_kmers import main as extract
import pandas as pd
import os
import sys
from tqdm import tqdm
from tabulate import tabulate


def split_fasta(file:os.PathLike,outfolder:os.PathLike)->None:

    with open(file) as handle:
        for values in SimpleFastaParser(handle):
            if not ".fa" in file: # just in case they don't add the extension
                out = open(outfolder+values[0]+".fa", 'w')
            else:
                out = open(outfolder+values[0], 'w')
            out.write('>'+values[0]+'\n'+values[1])
            out.close()

    return

if __name__ == "__main__":
    file = sys.argv[1]
    outfolder = sys.argv[2]
    split_fasta(file=file,outfolder=outfolder)