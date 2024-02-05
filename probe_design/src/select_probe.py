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
import shutil



def select_probe(currentfolder:os.PathLike = './data/',
                cutoff_cost:float = 1e5,      # max total cost of a probe. Exclude probe if even one oligo has a prohibitive cost
                cutoff_d_pc:int = 10,        # max distance between 2 consecutive oligos, as a % of the total probe length
                cutoff_d:int = 500           # max distance between 2 consecutive oligos, in nucleotides
                )->None:
    
    # independent variables moved to function arguments


    # retrieve probe queries
    selectedfolder = currentfolder + 'selected_probes/'
    try:
        os.mkdir(selectedfolder)
    except FileExistsError:
        pass
    # identify probe files
    pattern = currentfolder+"probe_candidates/**/probe_*.tsv"   #probelet
    #pattern = "data/**/probe*/**oligos.tsv"    #ifpd2 query
    #pattern = "data/**/*.best_probe.tsv"   #jupyternb
    filenames = glob.glob(pattern,recursive=True)

    roilist = currentfolder+'rois/all_regions.tsv'
    rdroi = pd.read_csv(roilist,sep="\t",header=0)


    probelist = pd.DataFrame(columns=['fullpath', 'folder', 'probe_set', 'roi', 'chr', 'probe_start', 'probe_end', \
            'nOligos', 'pw', 'probe_size', 'region_size', 'coverage', \
            'centrality', 'd_max_pcregion', 'd_max_pcprobe',\
            'd_mean', 'd_min', 'd_max', 'd_std', \
            'tm_range', 'tm_mean', 'tm_std', \
            'gc_range', 'gc_mean', 'gc_std', \
            'mean_closestMM', 'min_closestMM', 'max_closestMM', 'std_closestMM', \
            'mean_cumulMM', 'min_cumulMM', 'max_cumulMM', 'std_cumulMM', \
            'mean_oligo_cost', 'min_oligo_cost', 'max_oligo_cost', 'std_oligo_cost' ,\
            'sum_inv_cost'])
            #'mean_sumMM', 'min_sumMM', 'max_sumMM', 'std_sumMM'   ])
        

    # read probe files
    for file in tqdm(filenames, "Processing probe candidates"):
        rd = pd.read_csv(file,sep="\t")
        last = len(probelist)+1
        
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

        probelist.loc[last] = [file, filesplit[3], filesplit[4], rd.name[1], rd.chromosome[1], probe_start, probe_end, \
            len(rd), float(filesplit[4][-8:-4]), probe_end-probe_start+1, int(roi_end)-int(roi_start)+1, (probe_end-probe_start+1)/(int(roi_end)-int(roi_start)+1), \
            min(probe_center/roi_center,2-(probe_center/roi_center)), 100*(rd.start-rd.end.shift()).max()/(int(roi_end)-int(roi_start)+1), 100*(rd.start-rd.end.shift()).max()/(int(probe_end)-int(probe_start)+1),\
            (rd.start-rd.end.shift()).mean(), (rd.start-rd.end.shift()).min(), (rd.start-rd.end.shift()).max(), (rd.start-rd.end.shift()).std(), \
            max(rd.Tm)-min(rd.Tm), stat.mean(rd.Tm), stat.stdev(rd.Tm), \
            max(rd.gc_content)-min(rd.gc_content), stat.mean(rd.gc_content), stat.stdev(rd.gc_content), \
            stat.mean(rd.off_target_no), min(rd.off_target_no), max(rd.off_target_no), stat.stdev(rd.off_target_no), \
            stat.mean(rd.off_target_sum), min(rd.off_target_sum), max(rd.off_target_sum), stat.stdev(rd.off_target_sum), \
            stat.mean(rd.oligo_cost), min(rd.oligo_cost), max(rd.oligo_cost), stat.stdev(rd.oligo_cost), \
            inv_cost.sum()]
        # stat.mean(rd.off_target_sum), min(rd.off_target_sum), max(rd.off_target_sum), stat.stdev(rd.off_target_sum)]  
    
    probelist.loc[:,'fullpath'] = probelist.loc[:,'fullpath'].astype(str)

    # one selected probe per ROI
    for roi in tqdm(pd.unique(probelist.roi),"Selecting optimal probes..."):
        roiprobelist = probelist[(probelist.max_oligo_cost < cutoff_cost) & (probelist.d_max_pcprobe < cutoff_d_pc) & (probelist.d_max < cutoff_d) & (probelist.roi == roi)]
        filteredroilist = roiprobelist[roiprobelist.nOligos == roiprobelist.nOligos.max()]
        selectedprobeindex = filteredroilist[['pw']].idxmin()
        shutil.copy2(filteredroilist.loc[selectedprobeindex].fullpath.tolist()[0],selectedfolder)

    return

if __name__ == "__main__":
    # select_probe.py (cuttoff_d)
    if (len(sys.argv)==2):
        cutoff_d = sys.argv[1]
    else:
        cutoff_d = 500 # max distance between 2 consecutive oligos, in nucleotides
    select_probe(cutoff_d = cutoff_d)