#!/usr/bin/python

import numpy as np
from tabulate import tabulate as tab
import os
from tqdm import tqdm
import pandas as pd
import sys


if __name__ == '__main__':
    #syntax: ./reform_hush.py DNA/RNA/-RNA 80 22

    L = int(sys.argv[2])     #full oligo length
    l = int(sys.argv[3])
    types = {'DNA' : 'Reference', 'RNA' : 'RevCompl', '-RNA' : 'Reference'}
    suffix = types[sys.argv[1]]


    currentfolder = './data'

    roilist = currentfolder+'/rois/all_regions.tsv'
    rd = pd.read_csv(roilist,sep="\t",header=0)
    ROIcount = len(rd)

    infolder = currentfolder+'/candidates/'
    out = currentfolder+'/HUSH_candidates/'

    try:
        os.mkdir(out)
    except FileExistsError:
        print("Saving to existing 'HUSH_candidates' directory.")

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
            splitseq[2*k2+1] = splitseq[2*k2+1] + ", " + "111" + str(hdist_grouped_min[k2]) + "987" + str(hdist_grouped_sum[k2]) + "\n"
        
        o = open(out+filename,'w')
        outseq = ''.join(splitseq)
        o.write(outseq)
        o.close()
        f.close()
