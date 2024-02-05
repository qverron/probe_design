#!/usr/bin/python3

from pickle import FALSE, TRUE
import numpy as np
from tabulate import tabulate as tab
import os
from tqdm import tqdm
import pandas as pd
import sys

types = {'DNA' : 'Reference', 'RNA' : 'RevCompl', '-RNA' : 'Reference'}

def reform_hush(nt_type:str='DNA',
                sub:bool=False,L:int=80,l:int=22,
                currentfolder:os.PathLike = './data')->None:

    suffix = types[nt_type] # see types dict above
    roilist = currentfolder+'/rois/all_regions.tsv'
    rd = pd.read_csv(roilist,sep="\t",header=0)
    ROIcount = len(rd)

    infolder = currentfolder+'/candidates/'
    out = currentfolder+'/HUSH_candidates/'

    try:
        os.mkdir(out)
    except FileExistsError:
        print("Saving to existing 'HUSH_candidates' directory.")

    if (sub==TRUE):
        # each oligo gets an off-target score based on the min number of mismatch in any of the sublengths, 
        # and the sum of mismatches count across all sublength oligos

        for k1 in tqdm(range(ROIcount), 'Processing all ROI'):
            filename = 'roi_'+str(k1+1)+'.GC35to85_'+suffix+'.fa'
            fasta = infolder+filename
            hush = fasta+'.'+str(l)+'mers.nh.L'+str(l)+'.mindist.uint8'
            hdist = np.fromfile(hush,'uint8')

            n = len(hdist)
            #print('HUSH list: '+str(n))
            hdist_grouped_min = hdist.reshape([n//(L-l+1), (L-l+1)]).min(axis=1)
            hdist_grouped_sum = hdist.reshape([n//(L-l+1), (L-l+1)]).sum(axis=1)
            #print('Grouped length: '+str(len(hdist_grouped_min)))
                    
            f = open(fasta,'r')
            full = f.read()
            splitseq = full.splitlines()

            for k2 in range(len(hdist_grouped_min)):
                #print('Current row: '+str(k2))
                splitseq[2*k2] = splitseq[2*k2] + "\n"
                # /!\ UGLY SOLUTION
                # encode the two mismatch parameters (min mismatch hit and sum mismatch count) in a format compatible with ifpd2 db make
                splitseq[2*k2+1] = splitseq[2*k2+1] + ", " + "111" + str(hdist_grouped_min[k2]) + "987" + str(hdist_grouped_sum[k2]) + "\n"
                    
            o = open(out+filename,'w')
            outseq = ''.join(splitseq)
            o.write(outseq)
            o.close()
            f.close()  

    else:
        # each oligo gets an off-target score based on the min number of mismatches as found by nHUSH

        for k1 in tqdm(range(ROIcount), 'Processing all ROI'):
            filename = 'roi_'+str(k1+1)+'.GC35to85_'+suffix+'.fa'
            fasta = infolder+filename
            hush = fasta+'.nh.L'+str(L)+'.mindist.uint8'
            hdist = np.fromfile(hush,'uint8')

            n = len(hdist)                    
            f = open(fasta,'r')
            full = f.read()
            splitseq = full.splitlines()

            for k2 in range(n):
                #print('Current row: '+str(k2))
                splitseq[2*k2] = splitseq[2*k2] + "\n"
                # /!\ UGLY SOLUTION
                # encode the two mismatch parameters (min mismatch hit and sum mismatch count) in a format compatible with ifpd2 db make
                splitseq[2*k2+1] = splitseq[2*k2+1] + ", " + "111" + str(hdist[k2]) + "98799"  + "\n"
                    
            o = open(out+filename,'w')
            outseq = ''.join(splitseq)
            o.write(outseq)
            o.close()
            f.close()
    return

if __name__ == "__main__":
    #syntax: ./reform_hush.py DNA/RNA/-RNA 80 (22)

    nt_type = sys.argv[1] # DNA/RNA/-RNA
    print(f'Number of arguments: '+str(len(sys.argv)))
    if(len(sys.argv) == 3):
        # no sublength was used
        sub = FALSE
    elif(len(sys.argv) == 4):
        sub = TRUE
    else:
         print(f'Incorrect number of arguments. Exiting...')
         exit(-1)

    L = int(sys.argv[2])     #full oligo length
    print(f'Length: '+str(L))
    if (sub==TRUE):
        l = int(sys.argv[3]) # sublength
        print(f'Sublength: '+str(l))
        reform_hush(nt_type=nt_type,sub=sub,L=L,l=l)
    else:
        print(f'No sublength used')
        reform_hush(nt_type=nt_type,sub=sub,L=L)
    