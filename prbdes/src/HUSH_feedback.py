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
from itertools import compress

# for each result file from HUSH
# read all rows and fetch HUSH score
# read corresponding probe candidate tsv file
# attribute HUSH score to tsv file rows

if __name__ == '__main__':
    currentfolder = './data/'
    selectedfolder = currentfolder + 'selected_probes/'
    # identify probe files
    hushpattern = currentfolder+"selected_probes/query_*.out"   # HUSH validation output
    hushnames = glob.glob(hushpattern,recursive=True)

    dbpattern = currentfolder+"db_tsv/db.*.tsv"   # ROI oligo databases
    dbnames = glob.glob(dbpattern,recursive=True)
    if(len(sys.argv)>2):
        cutoff = sys.argv[1]
        wipe = sys.argv[2]
    if(len(sys.argv)>1):
        cutoff = sys.argv[1]
    else:    
        cutoff = 99         # max HUSH score

    # for each result file from HUSH
    # read all rows and fetch HUSH score
    for hushfile in hushnames:
        #print(f'Current: '+hushfile)
        # identify DB file
        filesplit = hushfile.split('/')
        basename = filesplit[len(filesplit)-1]
        basename = basename[12:basename.find('.')]      # extracted roi name
        dbtest = [db.find(basename)>0 for db in dbnames]
        roiDb = list(compress(dbnames,dbtest))[0]

        print(f'Currently processing '+basename+'.')

        hush_fasta = open(hushfile, 'r')  
        hushscores = []
        lines = hush_fasta.readlines()
        for line in lines:
            linesplit = line.split(sep=", ")
            if len(linesplit) > 1:
                hushscores.append(int(linesplit[1][:-1]))    

        filename = hushfile.replace('query_','')
        tsvfile = filename[:filename.find('.fa')]+".tsv"

        #print(f'TSV file: '+tsvfile)

        # retrieve oligos with poor HUSH score
        candidate = pd.read_csv(tsvfile,sep="\t",header=0)
        candidate = candidate.set_index('start')
        candidate['HUSH_score'] = hushscores
        exclude = candidate[candidate.HUSH_score > cutoff]

        finished = True

        if(len(exclude)>0):
            print(f'Excluding '+str(len(exclude))+' oligos from the database.')
            # attribute prohibitive escafish score in the ROI oligo database
            roioligos = pd.read_csv(roiDb,sep="\t",header=0)
            cols = list(roioligos.columns.values)       # save the order of the columns to restore it later
            roioligos = roioligos.set_index('start')
            roioligos.loc[exclude.index,'oligo_cost'] = 1e10

            # export the updated database in place
            roioligos = roioligos.reset_index()
            roioligos = roioligos[cols]
            roioligos.to_csv(roiDb,index=False,sep="\t")

            finished = False
        else:
            print(f'No oligos were excluded.')    

        if (not finished):
            print("Not completed. Run again!")

        
        
        
