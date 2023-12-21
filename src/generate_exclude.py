#!/usr/bin/python3
#set -e

# Generate default masked-out regions, corresponding to the region coordinates as provided by the user.
# Used to avoid excluding region-specific repetitive sequences when running (n)HUSH.

# No arguments required.

from cmath import exp
import numpy as np
from tabulate import tabulate as tab
import os
from tqdm import tqdm
import pandas as pd 
import sys
import contextlib
from tqdm import tqdm


if __name__ == '__main__':
    currentfolder = './data'

    roilist = currentfolder+'/rois/all_regions.tsv'
    rd = pd.read_csv(roilist,sep="\t",header=0)
    ROIcount = len(rd)

    out = currentfolder+'/exclude/'

    try:
        os.mkdir(out)
    except FileExistsError:
        print("The data/exclude folder already exists, risk of overwriting manual input! Exiting.")
    else:    
        
        # for all ROI in the region list, create an exclusion file with the region coordinates
        for roi in tqdm(range(ROIcount),'Generating exclusion maps...'):
            fixedchr = rd.ref[roi]+".chromosome."+rd.chrom[roi][3:]
            exporttable = pd.DataFrame(columns=['chrom','chromStart','chromEnd'])
            exporttable.loc[roi]=[fixedchr,rd.Window_start[roi],rd.Window_end[roi]]
            filename = 'excl_roi_'+str(rd.window_id[roi])+'.bed'
            exporttable.to_csv(out+filename,index=False,sep="\t") 
