#!/usr/bin/python3

from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd
import os
import glob
import sys
from tqdm import tqdm 
from tabulate import tabulate
import statistics as stat



currentfolder = './data/'
# identify probe files
pattern = currentfolder+"final_probes/probe_*.tsv"   #probelet
#pattern = "data/**/probe*/**oligos.tsv"    #ifpd2 query
#pattern = "data/**/*.best_probe.tsv"   #jupyternb
filenames = glob.glob(pattern,recursive=True)

roilist = currentfolder+'rois/all_regions.tsv'
rdroi = pd.read_csv(roilist,sep="\t",header=0)


table = pd.DataFrame(columns=['folder', 'probe_set', 'roi', 'chr', 'probe_start', 'probe_end', \
        'nOligos', 'pw', 'probe_size', 'region_size', 'coverage', \
        'centrality', 'd_max(%RegionSize)', \
        'd_mean', 'd_min', 'd_max', 'd_std', \
        'tm_range', 'tm_mean', 'tm_std', \
        'gc_range', 'gc_mean', 'gc_std', \
        'mean_consecOffTarget', 'min_consecOffTarget', 'max_consecOffTarget', 'std_consecOffTarget', \
        'mean_cumulMM', 'min_cumulMM', 'max_cumulMM', 'std_cumulMM', \
        'mean_oligo_cost', 'min_oligo_cost', 'max_oligo_cost', 'std_oligo_cost' ,\
        'sum_inv_cost'])
        #'mean_sumMM', 'min_sumMM', 'max_sumMM', 'std_sumMM'   ])
    

# read probe files
for file in tqdm(filenames, "Processing final probes"):
    rd = pd.read_csv(file,sep="\t")
    last = len(table)+1
    
    # find coordinates of the corresponding roi
    current_roi = rd.name[1][4:]
    for k in range(len(rdroi)):
        if str(rdroi.window_id[k]) == current_roi:
            current = k

    roi_start = rdroi.Window_start[current]
    roi_end = rdroi.Window_end[current]
    roi_center = ((roi_start+roi_end)/2)-roi_start  #adjusted roi center

# calculate summary measures and add as last row

    inv_cost = np.divide(np.ones(len(rd.name.tolist())),rd.oligo_cost)

# replace off_target_no by min_mismatch for new queries
    probe_start = min(rd.start[1:])
    probe_end = max(rd.end[1:])
    probe_center = ((probe_start+probe_end)/2)-roi_start    #adjusted probe center

    filesplit = file.split("/")

    table.loc[last] = [filesplit[2], filesplit[3], int(current_roi), rd.chromosome[1], probe_start, probe_end, \
        len(rd), filesplit[3][-8:-4], probe_end-probe_start+1, int(roi_end)-int(roi_start)+1, (probe_end-probe_start+1)/(int(roi_end)-int(roi_start)+1), \
        min(probe_center/roi_center,2-(probe_center/roi_center)), 100*(rd.start-rd.end.shift()).max()/(int(roi_end)-int(roi_start)+1), \
        (rd.start-rd.end.shift()).mean(), (rd.start-rd.end.shift()).min(), (rd.start-rd.end.shift()).max(), (rd.start-rd.end.shift()).std(), \
        max(rd.Tm)-min(rd.Tm), stat.mean(rd.Tm), stat.stdev(rd.Tm), \
        max(rd.gc_content)-min(rd.gc_content), stat.mean(rd.gc_content), stat.stdev(rd.gc_content), \
        stat.mean(rd.off_target_no), min(rd.off_target_no), max(rd.off_target_no), stat.stdev(rd.off_target_no), \
        stat.mean(rd.off_target_sum), min(rd.off_target_sum), max(rd.off_target_sum), stat.stdev(rd.off_target_sum), \
        stat.mean(rd.oligo_cost), min(rd.oligo_cost), max(rd.oligo_cost), stat.stdev(rd.oligo_cost), \
        inv_cost.sum()]
       # stat.mean(rd.off_target_sum), min(rd.off_target_sum), max(rd.off_target_sum), stat.stdev(rd.off_target_sum)]  
 

# append to existing probe summary file
output = currentfolder+"final_probes_summary.tsv"
sumExists = os.path.isfile(output)
#if (sumExists):
#    outtable = pd.read_csv(output,sep="\t")
#    newtable = pd.concat([outtable, table])
#else:
newtable = table   

newtable.sort_values(by=['roi'], inplace=True)

# write into probe summary file
newtable.to_csv(output,index=False,sep="\t")
