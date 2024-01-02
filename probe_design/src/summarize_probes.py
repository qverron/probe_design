from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd
import os
import glob
import sys
from tqdm import tqdm 
from tabulate import tabulate
import statistics as stat


def summarize_probes(currentfolder:os.PathLike = './data/')->None:
    # identify probe files
    pattern = currentfolder+"probe_candidates/**/probe_*.tsv"   #probelet
    #pattern = "data/**/probe*/**oligos.tsv"    #ifpd2 query
    #pattern = "data/**/*.best_probe.tsv"   #jupyternb
    filenames = glob.glob(pattern,recursive=True)

    roilist = currentfolder+'/rois/all_regions.tsv'
    rdroi = pd.read_csv(roilist,sep="\t",header=0)


    table = pd.DataFrame(columns=['folder', 'probe_set', 'roi', 'chr', 'start', 'end', \
            'nOligos', 'pw', 'size', 'coverage', \
            'centrality', \
            'd_mean', 'd_min', 'd_max', 'd_std', \
            'tm_range', 'tm_mean', 'tm_std', \
            'gc_range', 'gc_mean', 'gc_std', \
            'mean_closestMM', 'min_closestMM', 'max_closestMM', 'std_closestMM', \
            'mean_oligo_cost', 'min_oligo_cost', 'max_oligo_cost', 'std_oligo_cost' ])
            #'mean_sumMM', 'min_sumMM', 'max_sumMM', 'std_sumMM'   ])
        

    # read probe files
    for file in tqdm(filenames, "Processing probe candidates"):
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

    # replace off_target_no by min_mismatch for new queries
        probe_start = min(rd.start[1:])
        probe_end = max(rd.end[1:])
        probe_center = ((probe_start+probe_end)/2)-roi_start    #adjusted probe center

        filesplit = file.split("/")

        table.loc[last] = [filesplit[2], filesplit[3], rd.name[1], rd.chromosome[1], probe_start, probe_end, \
            len(rd), filesplit[2][-4:], probe_end-probe_start+1, (probe_end-probe_start+1)/(int(roi_end)-int(roi_start)), \
            min(probe_center/roi_center,2-(probe_center/roi_center)), \
            (rd.start-rd.end.shift()).mean(), (rd.start-rd.end.shift()).min(), (rd.start-rd.end.shift()).max(), (rd.start-rd.end.shift()).std(), \
            max(rd.Tm)-min(rd.Tm), stat.mean(rd.Tm), stat.stdev(rd.Tm), \
            max(rd.gc_content)-min(rd.gc_content), stat.mean(rd.gc_content), stat.stdev(rd.gc_content), \
            stat.mean(rd.off_target_no), min(rd.off_target_no), max(rd.off_target_no), stat.stdev(rd.off_target_no), \
            stat.mean(rd.oligo_cost), min(rd.oligo_cost), max(rd.oligo_cost), stat.stdev(rd.oligo_cost)]
        # stat.mean(rd.off_target_sum), min(rd.off_target_sum), max(rd.off_target_sum), stat.stdev(rd.off_target_sum)]  
    

    # append to existing probe summary file
    output = currentfolder+"probe_summary_probelet.tsv"
    sumExists = os.path.isfile(output)
    #if (sumExists):
    #    outtable = pd.read_csv(output,sep="\t")
    #    newtable = pd.concat([outtable, table])
    #else:
    newtable = table   

    newtable.sort_values(by=['probe_set'], inplace=True)

    # write into probe summary file
    newtable.to_csv(output,index=False,sep="\t")

    return

if __name__ == '__main__':
    summarize_probes()