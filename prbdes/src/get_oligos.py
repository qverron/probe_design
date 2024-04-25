#!/usr/bin/python3


# TO DO: 
# Switch to click arguments


from Bio.SeqIO.FastaIO import SimpleFastaParser
from ifpd2q.scripts.extract_kmers import main as extract
import pandas as pd
import os
import sys
from tqdm import tqdm 
from tabulate import tabulate

# Retrieve complete sequences in ROI
if __name__ == '__main__':
    #syntax: ./get_oligos.py DNA/RNA gcfilter extfolder
    gcfilter = 1
    extfolder = './data/'

    if(len(sys.argv) == 1):
        # no probe type was specified
        type = 'DNA'
    elif(len(sys.argv) == 2):
        type = sys.argv[1]
    elif(len(sys.argv) == 3):
        type = sys.argv[1]
        gcfilter = sys.argv[2]
    elif(len(sys.argv) == 4):
        type = sys.argv[1]
        gcfilter = sys.argv[2]
        extfolder = sys.argv[3]

    else:
         print(f'Incorrect number of arguments. Exiting...')
         exit(-1)


    if type=='DNA':
        roifile = extfolder+'rois/all_regions.tsv'
        ref = extfolder+'ref/'
        outseq = extfolder+'regions/'

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
        outcan = extfolder+'candidates/'

        try:
            os.mkdir(outcan)
        except FileExistsError:
            print("Saving to existing 'candidates' directory.")
            
        for k in range(len(rd)):
            fullseq = outseq+'roi_'+str(rd.window_id[k])+'.fa'
            if not os.path.isfile(fullseq):
                print('The FASTA sequence for ROI '+str(rd.window_id[k])+' is missing.')
                continue
            extract(fullseq,outcan,rd.length[k],gcfilter)

    elif type == 'RNA':
        roifile = extfolder+'rois/all_regions.tsv'
        ref = extfolder+'ref/'
        outseq = extfolder+'regions/'

        # ASSUME THAT THE TRANSCRIPT SEQUENCES HAVE ALREADY BEEN IDENTIFIED
        # implement direct transcript retrieval? 

        f = open(roifile)
        rd = pd.read_csv(f,sep="\t",header=0)              
                        
        # Divide into k-mers
        outcan = extfolder+'candidates/'

        try:
            os.mkdir(outcan)
        except FileExistsError:
            print("Saving to existing 'candidates' directory.")
            
        for k in range(len(rd)):
            fullseq = outseq+'roi_'+str(rd.window_id[k])+'.fa'
            if not os.path.isfile(fullseq):
                print('The FASTA sequence for ROI '+str(rd.window_id[k])+' is missing.')
                continue
            extract(fullseq,outcan,rd.length[k],gcfilter)     
