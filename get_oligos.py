#!/usr/bin/python

from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import os
import sys
from tqdm import tqdm 
from tabulate import tabulate

# Retrieve complete sequences in ROI

roifile = './data/rois/all_regions.tsv'
ref = './data/ref/'
outseq = './data/regions/'

try:
    os.mkdir(outseq)
except FileExistsError:
    print("Saving to existing 'regions' directory.")
    
f = open(roifile)
rd = pd.read_csv(f,sep="\t",header=0)

for k in tqdm(range(len(rd)),desc='Retrieving sequences for all ROIs'):
    if not pd.isnull(rd.ref[k]):
        chrnr = rd.chrom[k][3:]
        reffile = ref+str(rd.at[k,'ref'])+'.chromosome.'+chrnr+'.fa'

        with open(reffile) as handle:
            for values in SimpleFastaParser(handle):
                fullseq = values[1]
                seq = fullseq[rd.Window_start[k]-1:rd.Window_end[k]]         # shift index by 1 to match ref genome

                # export sequences
                out = open(outseq+'roi_'+str(rd.window_id[k])+'.fa','w')
                out.write('>ROI_'+str(rd.window_id[k])+' pos='+rd.chrom[k]+':'+str(rd.Window_start[k])+'-'+str(rd.Window_end[k])+'\n'+seq)
                out.close()       
                
                
# Divide into k-mers

from ifpd2_new.ifpd2.scripts.extract_kmers import main as extract

outcan = './data/candidates/'

try:
    os.mkdir(outcan)
except FileExistsError:
    print("Saving to existing 'candidates' directory.")
    
for k in range(len(rd)):
    fullseq = outseq+'roi_'+str(rd.window_id[k])+'.fa'
    if not os.path.isfile(fullseq):
        print('The FASTA sequence for ROI '+str(rd.window_id[k])+' is missing.')
        continue
    extract(fullseq,outcan,rd.length[k])