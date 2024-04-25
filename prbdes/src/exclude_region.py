#!/usr/bin/python3

# Excludes a certain region from a reference chromosome.
# Used to run HUSH on particularly repetitive sequences
# check for matching oligos outside of the region of interest and its repeats

# No arguments
# Exclude regions according to bedfile in the "./data/exclude" folder, called excl_roi_#
# Fetches the chromosome in the "./data/ref/" folder and generates a genome ref file with that region masked

from Bio.SeqIO.FastaIO import SimpleFastaParser
from ifpd2q.scripts.extract_kmers import main as extract
import pandas as pd
import os
import sys
from tqdm import tqdm
from tabulate import tabulate


if __name__ == '__main__':

    roifile = './data/rois/all_regions.tsv'     # proceed all ROIs as provided in region list
    ref = './data/ref/'
    excl = './data/exclude/'

    f = open(roifile)
    rd = pd.read_csv(f, sep="\t", header=0)

    for k in tqdm(range(len(rd)), desc='Processing all ROIs'):
        # fetch bed file with sequences to exclude
        bed = excl + 'excl_roi_' + str(rd.window_id[k]) + '.bed'

        if os.path.isfile(bed):   
            # load bed file corresponding to the current roi  
            bd = pd.read_csv(bed, delimiter='\t',header=0)

            # create a dictionary of all chromosomes
            dictseq = {}
            files = os.listdir(ref)
            for i in range(len(files)):
                with open(ref+files[i]) as handle:
                    if files[i][0:3] != "gen":
                        seq = ""
                        for values in SimpleFastaParser(handle):
                            seq = seq + values[1]  
                        dictseq[files[i]] = seq 
                    else:   # exclude genome ref files from dictionary
                        dictseq[files[i]] = ""    


            # edit chromosomes according to bed file
            for j in range(len(bd)):
                # in case the bed file refers to chromosomes not found in the ref folder
                if os.path.isfile(ref+bd.chrom[j]+".fa"):
                    # instead of cutting out the repeats, replace with NNNN (avoids index shifting issues)
                    current = dictseq[bd.chrom[j]+".fa"]
                    insert = "N"*(bd.chromEnd[j]-bd.chromStart[j]+1)
                    new = current[:bd.chromStart[j]] + insert + current[bd.chromEnd[j]+1:]      # indexing: T[a:b] a included b not
                    # replace in chromosome dictionary
                    dictseq[bd.chrom[j]+".fa"] = new
                else:
                    print("The chromosome "+bd.chrom[j]+".fa was not found in the ref folder (ROI "+str(rd.window_id[k])+"). Skipping.")    

            # concatenate chromosomes into new ref genome
            genome = ""
            for c in range(len(files)):
                genome = genome + dictseq[files[c]]

            # export new ref genome
            out = open(ref+'genome_roi_' + str(rd.window_id[k])+'.fa', 'w')
            out.write('>Genome_ROI_'+str(rd.window_id[k])+'\n'+genome)
            out.close()
