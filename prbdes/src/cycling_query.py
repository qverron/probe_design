#!/usr/bin/env python

from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd
import os
import glob
import sys
import click
import re
from tqdm import tqdm 
from tabulate import tabulate
import statistics as stat
import shutil
from itertools import compress
import subprocess
import logging
from datetime import datetime
import multiprocessing as mp
from joblib import Parallel, delayed
import joblib
import contextlib

pd.options.mode.chained_assignment = None  # default='warn'. Suppress SettingWithCopyWarning

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

@click.command(
    name="cycling_query",
    help="Generate an optimized probe set for each ROI."
)
@click.option('-s', '--strand', type=click.STRING)
@click.option('-L', '--length', type=click.INT)
@click.option('-m', '--mismatch', type=click.INT)
@click.option('-c', '--cutoff', type=click.INT)
@click.option('-t', '--threads', type=click.INT)
@click.option('-g', '--gap', type=click.INT)
@click.option('-gpercent', '--gappercent', type=click.INT)
@click.option('-stepdown', type=click.INT)
@click.option('-probe', type=click.INT)
@click.option('-start', type=click.INT)
@click.option('-end', type=click.INT)
@click.option('-step', type=click.INT)
@click.option('-greedy', is_flag=True)
@click.option('-excl', is_flag=True)
@click.option('-noquerylog', is_flag=True)

def output(strand, length, mismatch, cutoff, threads, gap, greedy,excl,noquerylog, gappercent = None, stepdown = None, probe=None,start=None, end=None, step=None):
        # initialization
    currentfolder = './data/'       # can be adapted so the code can be run in other folders
    logdir = currentfolder + "logfiles/"
    try:
        os.mkdir(logdir)
    except FileExistsError:
        pass

    now = datetime.now()
    nowstring = now.strftime("%Y_%m_%d_%H:%M:%S")
    logpath = logdir+'cycling_query_'+nowstring+'.log'

    logging.basicConfig(filename=logpath, level=logging.DEBUG)

    logging.info(f"Starting probe query: "+nowstring)

    logging.info(f"User-set parameters:")
    logging.info(f"FISH type                : {strand}")
    logging.info(f"Oligo length             : {length}")
    logging.info(f"HUSH mismatches          : {mismatch}")
    logging.info(f"Max nb of off-targets    : {cutoff}")
    if(excl):
        logging.info(f"Masking probe region from HUSH runs.")
    

    outprobes = currentfolder +"final_probes/"
    try:
        os.mkdir(outprobes)
    except FileExistsError:
        pass
    # TO DO: check already completed probes and skip

    # parameters
    cutoff_cost = 1e6       # max total cost of a probe. Exclude probe if even one oligo has a prohibitive cost
    if (not gappercent):
        cutoff_d_pc = 10
    else:
        cutoff_d_pc = gappercent    
    #cutoff_d_pc = 10        # max distance between 2 consecutive oligos, as a % of the total probe length
    cutoff_d = gap          # max distance between 2 consecutive oligos, in nucleotides 
    #cutoff_d = 500          # max distance between 2 consecutive oligos, in nucleotides 
    cutoff_oligo = 10       # max allowed cost for a single oligo
    if (not stepdown):
        stepdown = 1                # number of oligos to decrease size of probe with, if no valid probe could be found with the current size       

    logging.info(f"Encoded parameters:")
    logging.info(f"Max probe cost                               : {cutoff_cost}")
    logging.info(f"Max dist. between oligos (% probe length)    : {cutoff_d_pc}")
    logging.info(f"Max dist. between oligos (nt)                : {cutoff_d}")
    logging.info(f"Max single oligo cost                        : {cutoff_oligo}")

    # start process with full ROI list
    roilist = currentfolder+'rois/all_regions.tsv'
    rdroi = pd.read_csv(roilist,sep="\t",header=0)
    
    finished = False
    count = 1

    completelyfailed = []

    if ((not start) or (not end) or (not step)):          # check if the user has provided all three parameters to sweep different oligo numbers
        sweep = False  
    else:
        sweep = True
        oligorange = range(start,end,step)
    if (not probe):
        toprocess = rdroi
    else:
        toprocess = rdroi[rdroi.window_id == probe]      
          
    while (not finished):

        logging.warning(f"Probe generation round "+str(count)+".")

        print(f"Probe generation round "+str(count)+".")

        # EMPTY SELECTED_PROBES FOLDER + MAKE SURE PROBE QUERY FOLDERS DON'T COLLIDE
        try:
            shutil.rmtree(currentfolder + "probe_candidates")
        except:
            pass
        try:
            shutil.rmtree(currentfolder + "selected_probes") 
        except:
            pass       

        os.mkdir(currentfolder+"probe_candidates/")
        os.mkdir(currentfolder+"selected_probes/")
        
        toprocessRoi = toprocess.window_id.to_list()
        toprocessOligos = toprocess.window.to_list()

        # update the probe databases with removed oligos
        #filterdatabase(currentfolder,cutoff_oligo,toprocessRoi)    

        filterdatabase_par(currentfolder,cutoff_oligo,threads,toprocessRoi) 

        if excl:
            flag = "-e"
        else:
            flag=""    

        for n in tqdm(range(len(toprocess)),"Generating probe candidates..."):
            # retrieve ROI number from ROI name
            roinumber = toprocessRoi[n]
            oligos = toprocessOligos[n]
            timestamp = datetime.now()
            ts_string = timestamp.strftime("%Y%m%d_%H%M%S")

            if oligos <= 0:         # no probe was found at any length
                completelyfailed.append(roinumber)
                logging.warning(f"No probe could be found for ROI "+str(roinumber)+". Proceeding with the other probes.")
                continue

            # if the user provided start/end/step, use as range of oligo numbers to design probes for the first time!
            if sweep:
                querylogpath = logdir + "query_roi_"+str(roinumber)+"_oli_sweep_round_"+str(count)+"_" + ts_string + ".txt"
                timeout=9999999
                with tqdm_joblib(tqdm(desc="Sweeping oligo counts in region "+str(roinumber), total=len(oligorange))) as progress_bar:
                    Parallel(n_jobs=threads,timeout=timeout)(delayed(probequery)(length,strand,roinumber,oligos,querylogpath,greedy,noquerylog) for oligos in oligorange)
            else:
                querylogpath = logdir + "query_roi_"+str(roinumber)+"_oli_"+str(oligos)+"_round_"+str(count)+"_" + ts_string + ".txt"
                # use as input for probe query (only process remaining ROIs)
                probequery(length,strand,roinumber,oligos,querylogpath,greedy,noquerylog)

        # select best probes
        print(f"Selecting probes...")
        if(sweep):
            selection = selectprobes(currentfolder, toprocessRoi, [start]*len(toprocessRoi), cutoff_cost, cutoff_d, cutoff_d_pc) 
        else:
            selection = selectprobes(currentfolder, toprocessRoi, toprocessOligos, cutoff_cost, cutoff_d, cutoff_d_pc) 
    
        failedlist = selection[selection.success == 0].index.to_list()
        logging.info(f'Length of failedlist: '+str(len(failedlist)))
        
        selectedlist = selection[selection.success == 1].index.to_list()

        if(len(selectedlist)>0):
            # check off-target homology with HUSH
            timestamp = datetime.now()
            ts_string = timestamp.strftime("%Y%m%d_%H%M%S")
            hushlogpath = logdir + "hush_roi_round_"+str(count)+"_" + ts_string + ".txt"
            print(f"Checking the oligos with (old)HUSH...")
            with open(hushlogpath,'w') as f:
#            subprocess.run("./validation_oldHUSH_BLAST.sh -L "+str(length)+" -m "+str(mismatch)+" -t "+str(threads)+flag+" > "+hushlogpath, shell=True,check=True)
                subprocess.run(["prbdes","validation_oldHUSH_BLAST", '-L', str(length),"-m",str(mismatch),"-t",str(threads),flag],stdout=f)
            print(f"Removing poor oligos from database")

        # apply results from HUSH to exclude poor oligos
        rerunlist = feedback(currentfolder,outprobes,count,cutoff,logpath)

        combinedlist = np.unique(failedlist+rerunlist)

        if (len(combinedlist) == 0):
            finished = True
            logging.info(f"Done! :)")
            if(len(completelyfailed)>0):
                logging.info(f"No probe could be found for the following regions: "+''.join(str(e)+", " for e in completelyfailed)+".")
            break
        else:
            print(f""+str(len(combinedlist))+" probes need to be re-run.")
            
            # set the correct oligo numbers for the next iteration
            #  
            # if a probe was selected but rejected by HUSH, re-use the selected number of oligos
            maskrerun = [(toprocess.loc[k,'window_id'] in rerunlist) for k in toprocess.index.to_list()] 
            rerunrois = toprocess[maskrerun]
            if (len(rerunrois)>0):
                toprocess["window"][maskrerun] = [int(selection.loc[k,'oligos']) for k in rerunrois.window_id.to_list()]

            # for probes for which no valid probe could be constructed, reduce number of oligos for next iteration
            maskfailed = [(toprocess.loc[k,'window_id'] in failedlist) for k in toprocess.index.to_list()] 
            failedrois = toprocess[maskfailed]
            if (len(failedrois)>0):
                toprocess["window"][maskfailed] = [int(selection.loc[k,'oligos']-stepdown) for k in failedrois.window_id.to_list()]

            # only keep rois that need to be re-run
            maskcombined = [(toprocess.loc[k,'window_id'] in combinedlist) for k in toprocess.index.to_list()] 
            toprocess = toprocess[maskcombined]

            count = count+1  
            sweep = False   


# -----------------------------------------------------------------------------------------------------------------------      
# -----------------------------------------------------------------------------------------------------------------------            


def selectprobes(input_folder, toprocessroi, toprocessoligos, cutoff_cost, cutoff_d, cutoff_d_pc):

    # retrieve probe queries
    selectedfolder = input_folder + 'selected_probes/'
    try:
        os.mkdir(selectedfolder)
    except FileExistsError:
        pass
    # identify probe files
    pattern = input_folder+"probe_candidates/**/probe_*.tsv"   #probelet
    #pattern = "data/**/probe*/**oligos.tsv"    #ifpd2 query
    #pattern = "data/**/*.best_probe.tsv"   #jupyternb
    filenames = glob.glob(pattern,recursive=True)

    roilist = input_folder+'rois/all_regions.tsv'
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
    if len(filenames)>0:
        for file in tqdm(filenames, "Summarizing probe candidates..."):
            rd = pd.read_csv(file,sep="\t")
            last = len(probelist)+1
            
            # find coordinates of the corresponding roi
            current_roi = rd.name[0][4:]
            for k in range(len(rdroi)):
                if str(rdroi.window_id[k]) == current_roi:
                    current = k

            roi_start = rdroi.Window_start[current]
            roi_end = rdroi.Window_end[current]
            roi_center = ((roi_start+roi_end)/2)-roi_start  #adjusted roi center

        # calculate summary measures and add as last row

            inv_cost = np.divide(np.ones(len(rd.name.tolist())),rd.oligo_cost)

        # replace off_target_no by min_mismatch for new queries
            probe_start = min(rd.start[:])
            probe_end = max(rd.end[:])
            probe_center = ((probe_start+probe_end)/2)-roi_start    #adjusted probe center

            filesplit = file.split("/")         # split numbers will have to be adjusted for running in different folders

            probelist.loc[last] = [file, filesplit[3], filesplit[4], int(rd.name[0][4:]), rd.chromosome[0], probe_start, probe_end, \
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

    selection = pd.DataFrame(columns=['roi','oligos','pw','success'])   # summary of selection results
    selection.set_index('roi',inplace=True)

    # one selected probe per ROI
    for n in tqdm(range(len(toprocessroi)),"Selecting optimal probes..."):
        roi = toprocessroi[n]
        minoligo = toprocessoligos[n]
        try:
            roiprobelist = probelist[(probelist.max_oligo_cost < cutoff_cost) & (probelist.d_max_pcprobe < cutoff_d_pc) & (probelist.d_max < cutoff_d) & (probelist.roi == roi)]
            filteredroilist = roiprobelist[roiprobelist.nOligos == roiprobelist.nOligos.max()]
            selectedprobeindex = filteredroilist[['pw']].idxmin()
            shutil.copy2(filteredroilist.loc[selectedprobeindex].fullpath.tolist()[0],selectedfolder)
            selection.loc[roi] = [roiprobelist.nOligos.max(),roiprobelist.loc[selectedprobeindex,'pw'].to_list()[0],1]
        except:
            logging.info(f"No valid probe was found for ROI #"+str(roi)+".")
            if (len(probelist[probelist.roi == roi]) > 0):
                selection.loc[roi] = [probelist[probelist.roi == roi].nOligos.min(),0,0]    # smallest successful number of oligos
            else:
                selection.loc[roi] = [minoligo,0,0]   # min oligo count that was unsuccessfully tested   

    return selection


# -----------------------------------------------------------------------------------------------------------------------      
# -----------------------------------------------------------------------------------------------------------------------            


def feedback(currentfolder,outfolder,count,cutoff,logpath):
    selectedfolder = currentfolder + 'selected_probes/'
    # identify probe files
    hushpattern = currentfolder+"selected_probes/query_*.out"   # HUSH validation output
    hushnames = glob.glob(hushpattern,recursive=True)

    logging.info(f'Length of hushnames: '+str(len(hushnames)))
    #logging.info(f'Hushnames: '+' '.join(str(l) for l in hushnames))  

    dbpattern = currentfolder+"db_tsv/db.*.tsv"   # ROI oligo databases
    dbnames = glob.glob(dbpattern,recursive=True)

    rerunlist = []

    # for each result file from HUSH
    # read all rows and fetch HUSH score
    if(len(hushnames)>0):
        for hushfile in hushnames:

            #logging.info(f'Current hush file: '+hushfile)

            # identify DB file
            filesplit = hushfile.split('/')
            basename = filesplit[len(filesplit)-1]
            roiname = basename[12:basename.find('.')]      # extracted roi name
            dbtest = [db.find("."+roiname+".")>0 for db in dbnames]
            roiDb = list(compress(dbnames,dbtest))[0]      # identify the corresponding database file

            logging.info(f'Currently processing '+roiname+', round #'+str(count)+".")

            hush_fasta = open(hushfile, 'r')  
            hushscores = []
            lines = hush_fasta.readlines()
            for line in lines:
                linesplit = line.split(sep=", ")
                if len(linesplit) > 1:
                    hushscores.append(int(linesplit[1][:-1]))    

            filename = hushfile.replace('query_','')
            tsvfile = filename[:filename.find('.fa')]+".tsv"

            # retrieve oligos with poor HUSH score
            candidate = pd.read_csv(tsvfile,sep="\t",header=0)
            candidate = candidate.set_index('start')
            candidate['HUSH_score'] = hushscores
            exclude = candidate[candidate.HUSH_score > cutoff]

            if(len(exclude)>0):
                logging.info(f'Excluding '+str(len(exclude))+' oligos from the database.')
                # attribute prohibitive escafish score in the ROI oligo database
                roioligos = pd.read_csv(roiDb,sep="\t",header=0)
                cols = list(roioligos.columns.values)       # save the order of the columns to restore it later
                roioligos = roioligos.set_index('start')
                roioligos.loc[exclude.index,'oligo_cost'] = 1e10

                # export the updated database in place
                roioligos = roioligos.reset_index()
                roioligos = roioligos[cols]
                roioligos.to_csv(roiDb,index=False,sep="\t")

                rerunlist.append(int(roiname[4:]))  # the probe will have to be queried again from the updated oligo database

            else:
                logging.info(f'No oligos were excluded.')  
                # move probe to final selection folder
                tsvname = basename[6:basename.find('.fa')]+".tsv"
                shutil.move(selectedfolder + tsvname,outfolder)
                shutil.move(selectedfolder + basename,outfolder)    # also keep the .out file (all other files will be deleted)      
                    

    logging.info(f'Length of rerunlist: '+str(len(rerunlist)))
    logging.info(f'Rerunlist: '+' '.join(str(l) for l in rerunlist))

    return rerunlist


# -----------------------------------------------------------------------------------------------------------------------      
# -----------------------------------------------------------------------------------------------------------------------            


def filterdatabase(currentfolder,cutoff,toprocessRoi):
    # generate a filtered copy of the oligo database for each ROI
    # filter by removing oligos over a certain threshold cost

    dbfolder = currentfolder + 'db_tsv/'
    # identify database files
    dbpattern = dbfolder+"db.*.tsv"   # HUSH validation output
    dbnames = glob.glob(dbpattern,recursive=True)

    dbparse = [re.split('roi_',db,maxsplit=1)[1] for db in dbnames]
    dbrois = [int(re.split('\.',db,maxsplit=1)[0]) for db in dbparse]

    filtereddblist = [dbnames[k] for k in range(len(dbrois)) if dbrois[k] in toprocessRoi]

    # for each database file
    # open the updated (full) database and export a filtered version without discarded oligos       

    for db in tqdm(filtereddblist,"Preparing filtered oligo databases..."):
        # retrieve oligos with poor HUSH score
        oligodb = pd.read_csv(db,sep="\t",header=0)
        filtereddb = oligodb[oligodb.oligo_cost < cutoff]
        filtereddb.to_csv(db+".filt",index=False,sep="\t")

def filterdatabase_par(currentfolder,cutoff,threads,toprocessRoi):
    # generate a filtered copy of the oligo database for each ROI
    # filter by removing oligos over a certain threshold cost

    dbfolder = currentfolder + 'db_tsv/'
    # identify database files
    dbpattern = dbfolder+"db.*.tsv"   # HUSH validation output
    dbnames = glob.glob(dbpattern,recursive=True)

    dbparse = [re.split('roi_',db,maxsplit=1)[1] for db in dbnames]
    dbrois = [int(re.split('\.',db,maxsplit=1)[0]) for db in dbparse]

    filtereddblist = [dbnames[k] for k in range(len(dbrois)) if dbrois[k] in toprocessRoi]

    # for each database file
    # open the updated (full) database and export a filtered version without discarded oligos
    with tqdm_joblib(tqdm(desc="Filtering oligo databases", total=len(filtereddblist))) as progress_bar:
        Parallel(n_jobs=threads)(delayed(filterdb)(db,cutoff) for db in filtereddblist)
 

def filterdb(db,cutoff):
     # retrieve oligos with poor HUSH score
    oligodb = pd.read_csv(db,sep="\t",header=0)
    filtereddb = oligodb[oligodb.oligo_cost < cutoff]
    filtereddb.to_csv(db+".filt",index=False,sep="\t")       

# -----------------------------------------------------------------------------------------------------------------------      
# -----------------------------------------------------------------------------------------------------------------------            


def probequery(length,strand,roi,oligos,logpath,greedy,noquerylog):
    suffix = ""
    if(greedy): 
        suffix = "-g"
    if noquerylog:    
        #subprocess.run("./probe-query.sh -s "+strand+" -e "+str(roi)+" -o "+str(oligos)+suffix+"> /dev/null 2>&1", shell=True)
        subprocess.run(["prbdes","probe-query","-s",strand,"-e",str(roi),"-o",str(oligos),suffix], stdout=None)
    else:
        with open(logpath,'w') as f:
        #subprocess.run("./probe-query.sh -s "+strand+" -e "+str(roi)+" -o "+str(oligos)+suffix+" > "+logpath+" 2>&1", shell=True)
            subprocess.run(["prbdes","probe-query","-s",strand,"-e",str(roi),"-o",str(oligos),suffix],stderr=subprocess.STDOUT,stdout=f)

# -----------------------------------------------------------------------------------------------------------------------      
# -----------------------------------------------------------------------------------------------------------------------            


if __name__ == '__main__':
    output()
    logging.shutdown()
