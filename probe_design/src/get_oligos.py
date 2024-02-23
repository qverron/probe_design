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
def get_oligos(nt_type:str='DNA',
               gcfilter:bool = 1, 
               extfolder:os.PathLike = './data/',
               reffile:os.PathLike|None= None,
               )->None:

    if nt_type=='DNA':
        roifile = os.path.join(extfolder,'rois/all_regions.tsv')
        ref = os.path.join(extfolder,'ref/')
        outseq = os.path.join(extfolder,'regions/')

        try:
            os.mkdir(outseq)
        except FileExistsError:
            print("Saving to existing 'regions' directory.")
            
        f = open(roifile)
        rd = pd.read_csv(f,sep="\t",header=0)

        for k in tqdm(range(len(rd)),desc='Retrieving sequences for all ROIs'):
            if not pd.isnull(rd.ref[k]):
                chrnr = rd.chrom[k][3:] # chromosome number

                if reffile is None: # if no reference file is specified, use the one in the rois file
                    reffile = os.path.join(ref,f"{rd.at[k,'ref']}.chromosome.{chrnr}.fa")

                with open(reffile) as handle:
                    for values in SimpleFastaParser(handle):
                        fullseq = values[1]
                        seq = fullseq[rd.Window_start[k]-1:rd.Window_end[k]]         # shift index by 1 to match ref genome

                        # export sequences
                        out = open(os.path.join(outseq,'roi_'+str(rd.window_id[k])+'.fa'),'w')
                        out.write('>ROI_'+str(rd.window_id[k])+' pos='+rd.chrom[k]+':'+str(rd.Window_start[k])+'-'+str(rd.Window_end[k])+'\n'+seq)
                        out.close()       
                        
                        
        # Divide into k-mers
        outcan = os.path.join(extfolder,'candidates/')

        try:
            os.mkdir(outcan)
        except FileExistsError:
            print("Saving to existing 'candidates' directory.")
            
        for k in range(len(rd)):
            fullseq = os.path.join(outseq,f'roi_{rd.window_id[k]}.fa')
            if not os.path.isfile(fullseq):
                print('The FASTA sequence for ROI '+str(rd.window_id[k])+' is missing.')
                continue
            extract(fullseq,outcan,rd.length[k],gcfilter)

    elif type == 'RNA':
        roifile = os.path.join(extfolder,'rois/all_regions.tsv')
        ref = os.path.join(extfolder,'ref/')
        outseq = os.path.join(extfolder,'regions/')

        # ASSUME THAT THE TRANSCRIPT SEQUENCES HAVE ALREADY BEEN IDENTIFIED
        # implement direct transcript retrieval? 

        f = open(roifile)
        rd = pd.read_csv(f,sep="\t",header=0)              
                        
        # Divide into k-mers
        outcan = os.path.join(extfolder,'candidates/')

        try:
            os.mkdir(outcan)
        except FileExistsError:
            print("Saving to existing 'candidates' directory.")
            
        for k in range(len(rd)):
            fullseq = os.path.join(outseq,f"roi_{rd.window_id[k]}.fa")
            if not os.path.isfile(fullseq):
                print(f'The FASTA sequence for ROI {rd.window_id[k]} is missing.')
                continue
            extract(fullseq,outcan,rd.length[k],gcfilter)
    else:
        print(f'Invalid nt_type: {nt_type}. Exiting...')
        exit(-1)

    return   

if __name__ == '__main__':
    #syntax: ./get_oligos.py DNA/RNA gcfilter extfolder

    if(len(sys.argv) == 1):
        # no probe type was specified
        get_oligos(nt_type='DNA') # default is nt_type = 'DNA'
    elif(len(sys.argv) == 2):
        get_oligos(nt_type = sys.argv[1]) # default is nt_type = 'DNA'
    elif(len(sys.argv) == 3):
        get_oligos(nt_type = sys.argv[1],
                   gcfilter= sys.argv[2]) # gcfilter = 1 (True) by default
    elif(len(sys.argv) == 4):
        get_oligos(nt_type = sys.argv[1],
                   gcfilter = sys.argv[2],
                   extfolder = sys.argv[3]) # extraction folder = './data/' by default

    else:
         print(f'Incorrect number of arguments. Exiting...')
         exit(-1)
    