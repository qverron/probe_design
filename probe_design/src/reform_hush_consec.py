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
                maxccmatch = 0, # largest possible match,
                currentccmatch = 0,  # current match block,
                prev = 1, # pointer to previous sublength
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



def reform_hush_consec(nt_type:str='DNA',
                       sub:bool=False,
                       L:int=80,
                       l:int=22,
                       currentfolder:os.PathLike = './data',
                       )->None:


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

    if (sub==TRUE):
        
        for roi in tqdm(range(ROIcount),'Processing all regions'):
            ## each oligo gets an off-target score based on the max number of consecutive perfect matches (0 in nHUSH results)
            # and the resulting sum of such consecutive blocks
            # one 0 => one sublength perfect match
            # each consecutive 0 => +1 to length
            # count each block separately and sum the results
            
            filename = 'roi_'+str(rd.window_id[roi])+'.GC35to85_'+suffix+'.fa'
            fasta = infolder+filename
            hush = fasta+'.'+str(l)+'mers.nh.L'+str(l)+'.mindist.uint8'
            hdist = np.fromfile(hush,'uint8')

            n = len(hdist)

            hdist_grouped = hdist.reshape([n//(L-l+1), (L-l+1)])
            

            with tqdm_joblib(tqdm(desc="Processing all oligos in region "+str(rd.window_id[roi]), total=len(hdist_grouped))) as progress_bar:
                results = Parallel(n_jobs=40)(delayed(consecblock)(oligo,L,l,hdist_grouped) for oligo in range(len(hdist_grouped)))
            #for oligo in range(590):
            #    consecblock(oligo,L,l,hdist_grouped)

            #time.sleep(20000)

            # write results into fasta file  
            f = open(fasta,'r')
            full = f.read()
            splitseq = full.splitlines()

            for k in range(len(hdist_grouped)):
                #print('Current row: '+str(k2))
                splitseq[2*k] = splitseq[2*k] + "\n"
                # /!\ UGLY SOLUTION
                # encode the two mismatch parameters (min mismatch hit and sum mismatch count) in a format compatible with ifpd2 db make
                splitseq[2*k+1] = splitseq[2*k+1] + ", " + str(results[k]) + "\n"
                    
            o = open(out+filename,'w')
            outseq = ''.join(splitseq)
            o.write(outseq)
            o.close()
            f.close()  

    # else:
    #     # each oligo gets an off-target score based on the min number of mismatches as found by nHUSH

    #     for k1 in tqdm(range(ROIcount), 'Processing all ROI'):
    #         filename = 'roi_'+str(k1+1)+'.GC35to85_'+suffix+'.fa'
    #         fasta = infolder+filename
    #         hush = fasta+'.nh.L'+str(L)+'.mindist.uint8'
    #         hdist = np.fromfile(hush,'uint8')

    #         n = len(hdist)                    
    #         f = open(fasta,'r')
    #         full = f.read()
    #         splitseq = full.splitlines()

    #         for k2 in range(n):
    #             #print('Current row: '+str(k2))
    #             splitseq[2*k2] = splitseq[2*k2] + "\n"
    #             # /!\ UGLY SOLUTION
    #             # encode the two mismatch parameters (min mismatch hit and sum mismatch count) in a format compatible with ifpd2 db make
    #             splitseq[2*k2+1] = splitseq[2*k2+1] + ", " + "111" + str(hdist[k2]) + "98799"  + "\n"
                    
    #         o = open(out+filename,'w')
    #         outseq = ''.join(splitseq)
    #         o.write(outseq)
    #         o.close()
    #         f.close()         
    return

if __name__ == "__main__":
    # for Command Line use
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
        l = int(sys.argv[3])
        print(f'Sublength: '+str(l))
        reform_hush_consec(nt_type=nt_type,sub=sub,L=L,l=l)
    else:
        print(f'No sublength used')
        reform_hush_consec(nt_type=nt_type,sub=sub,L=L)
    
    


    