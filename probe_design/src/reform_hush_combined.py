#!/usr/bin/python3

from pickle import FALSE, TRUE
import numpy as np
import multiprocessing as mp
from joblib import Parallel, delayed
from tabulate import tabulate as tab
import os
from tqdm import tqdm
import pandas as pd 
import sys
import contextlib
import joblib
import time
from tqdm import tqdm

types = {'DNA' : 'Reference', 'RNA' : 'RevCompl', '-RNA' : 'Reference'}

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()



def consecblock(oligo,L,l,hdist_grouped,
                maxccmatch:int=0, # largest possible match
                currentccmatch:int=0, # current match block
                prev:int=1, # pointer to previous sublength
                )->int:
    # find consecutive 0 in each L-mer

    for sub in range(L-l+1):         # check the homology of all consecutive l-mers
        current = hdist_grouped[oligo,sub]
        if(current == 0 and prev>0):     # first zero in a new block
            currentccmatch = l 
        elif(current == 0 and prev==0):     # match block continues
            currentccmatch = currentccmatch + 1
        elif(current > 0 and prev==0):       # end of previous block    
            maxccmatch = max([maxccmatch,currentccmatch])
            currentccmatch = 0
        else:
            currentccmatch = 0

        #print(f'Oligo: '+str(oligo)+', current HUSH: '+str(current)+', current block: '+str(currentccmatch)+', max: '+str(maxccmatch)+', prev: '+str(prev))
        prev = current
    
    maxccmatch = max([maxccmatch,currentccmatch])
    
    return maxccmatch


def reform_hush_combined(nt_type:str='DNA',
                         L:int=40,
                         l:int=22,
                         currentfolder:os.PathLike = './data',
                         until:int=3)->None:
    
    suffix = types[nt_type]
    roilist = currentfolder+'/rois/all_regions.tsv'
    rd = pd.read_csv(roilist,sep="\t",header=0)
    ROIcount = len(rd)

    infolder = currentfolder+'/candidates/'
    out = currentfolder+'/HUSH_candidates/'

    try:
        os.mkdir(out)
    except FileExistsError:
        print("Saving to existing 'HUSH_candidates' directory.")

    ## each oligo gets an off-target score based on:
    # 1. The longest consecutive perfect off-target match (consecutive sublength oligos returning 0 in nHUSH results)
    # 2. The sum of mismatch counts in sublength oligos. Challenging to interpret but higher for central mismatches than closer to the edges.
    # Encode both scores to retrieve them separately after building db.

    for roi in tqdm(range(ROIcount),'Processing all regions'):
        
        filename = 'roi_'+str(rd.window_id[roi])+'.GC35to85_'+suffix+'.fa'
        fasta = infolder+filename
        hush = fasta+'.'+str(l)+'mers.nh.L'+str(l)+'.mindist.uint8'
        hdist = np.fromfile(hush,'uint8')

        # correct for aberrant values after nHUSH.
        hdist[hdist>until] = until+1
        n = len(hdist)

        hdist_grouped = hdist.reshape([n//(L-l+1), (L-l+1)])
        hdist_grouped_sum = hdist.reshape([n//(L-l+1), (L-l+1)]).sum(axis=1)            

        with tqdm_joblib(tqdm(desc="Processing all oligos in region "+str(rd.window_id[roi]), total=len(hdist_grouped))) as progress_bar:
            results = Parallel(n_jobs=40)(delayed(consecblock)(oligo,L,l,hdist_grouped) for oligo in range(len(hdist_grouped)))

        # write results into fasta file  
        f = open(fasta,'r')
        full = f.read()
        splitseq = full.splitlines()

        for k in range(len(hdist_grouped)):
            splitseq[2*k] = splitseq[2*k] + "\n"
            # /!\ UGLY SOLUTION
            # encode the two mismatch parameters (min mismatch hit and sum mismatch count) in a format compatible with ifpd2 db make
            splitseq[2*k+1] = splitseq[2*k+1] + ", " + "111" + str(results[k]) + "987" + str(hdist_grouped_sum[k]) + "\n"
                
        o = open(out+filename,'w')
        outseq = ''.join(splitseq)
        o.write(outseq)
        o.close()
        f.close()

        return

if __name__ == "__main__":
    #syntax: ./reform_hush_combined.py DNA/RNA/-RNA 40 22 3 

    print(f'Number of arguments: '+str(len(sys.argv)))
    if(len(sys.argv) < 5):
        print(f'The combined scoring approach requires nHUSH to be run using sublength oligos.\n')
        print(f'Required arguments: [DNA/RNA/-RNA] [length] [sublength] [until].\n')
        print(f'With "until" the same parameter as used when running nHUSH.\n')
        print(f'Incorrect number of arguments. Exiting...')
        exit(-1)

    # types = {'DNA' : 'Reference', 'RNA' : 'RevCompl', '-RNA' : 'Reference'}
    # suffix = types[sys.argv[1]]

    L = int(sys.argv[2])     #full oligo length
    print(f'Length: '+str(L))
    l = int(sys.argv[3])    # sublength used for nHUSH
    print(f'Sublength: '+str(l))
    until = int(sys.argv[4])
    print(f'HUSH was run until '+str(until)+' mismatches.')
    # call the function here
    reform_hush_combined(nt_type=sys.argv[1],L=L,l=l,until=until)  
