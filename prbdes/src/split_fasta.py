#!/usr/bin/python3

from Bio.SeqIO.FastaIO import SimpleFastaParser
from ifpd2q.scripts.extract_kmers import main as extract
import pandas as pd
import os
import sys
from tqdm import tqdm
from tabulate import tabulate


if __name__ == '__main__':
    file = sys.argv[1]
    outfolder = sys.argv[2]

    with open(file) as handle:
        for values in SimpleFastaParser(handle):
            out = open(outfolder+values[0]+".fa", 'w')
            out.write('>'+values[0]+'\n'+values[1])
            out.close()

