#!/usr/bin/python3

# GG's code sniplet (R)
#max_dg_fraction = 0.5
#max_off_targets = 99

#db[, off_target_score := 1-off_target_no/max_off_targets]
#db[off_target_score <= 0, off_target_score := 0]
#
#db[ss_dG/Tm_dG > max_dg_fraction, ss_score := 0]
#db[ss_dG/Tm_dG <= max_dg_fraction, ss_score := 1-ss_dG/Tm_dG/max_dg_fraction]
#db[ss_dG >= 0, ss_score := 1]
#
#db[, score := (off_target_score+ss_score)/2]

import sys
import re
from tqdm import tqdm 

def score_gg(d):
    # Extract variables from the dictionary
    ss_dG = float(d['ss_dG'])
    off_target_no = int(d['off_target_no'])
    Tm_dG = float(d['Tm_dG'])
    # Parameters
    max_dg_fraction = 0.5
    max_off_targets = 99
    # Calculate the score
    off_target_score = 1-off_target_no/max_off_targets
    if off_target_score < 0:
        off_target_score = 0
    ss_score = 0
    if(ss_dG/Tm_dG <= max_dg_fraction):
        ss_score = 1-ss_dG/Tm_dG/max_dg_fraction

    if ss_dG >= 0:
        ss_score = 1
    score = 1-(off_target_score + ss_score)/2
    return score

def score_gg_nhush(d):
    # Extract variables from the dictionary
    ss_dG = float(d['ss_dG'])
    off_target_no = int(d['off_target_no'])
    Tm_dG = float(d['Tm_dG'])
    # Parameters
    max_dg_fraction = 0.5
    max_off_targets = 99
    # Calculate the score
    # 0 -> 1, 1->1/2, 2->1/3 etc
    off_target_score = 1/(1+off_target_no)
    if off_target_score < 0:
        off_target_score = 0
    ss_score = 0
    if(ss_dG/Tm_dG <= max_dg_fraction):
        ss_score = 1-ss_dG/Tm_dG/max_dg_fraction

    if ss_dG >= 0:
        ss_score = 1
    #score = (off_target_score + ss_score)/2
    score = 1-0.5*ss_score + 0.5*off_target_score

    return score

def score_q_nhush(d):
    # Extract variables from the dictionary
    ss_dG = float(d['ss_dG'])
    off_target_no = int(d['off_target_no'])
    off_target_sum = int(d['off_target_sum'])
    Tm_dG = float(d['Tm_dG'])
    Tm = float(d['Tm'])
    seq = str(d['sequence'])
    # Parameters
    max_dg_fraction = 0.3
    max_offtarget = 4   # correction for aberrant values after mindist
    # Calculate the score
    # 0 -> 1, 1->1/2, 2->1/3 etc
    if(off_target_no>max_offtarget):
        off_target_score = 1e99
    else:    
        off_target_score = 1/(1+off_target_no)**2
    if(off_target_sum==0):
        sum_target_score = 1e99
    else:
        sum_target_score = 10/(1+off_target_sum)**(1/2)      
    ss_score = 1
    if(ss_dG/Tm_dG <= max_dg_fraction):
        ss_score = ss_dG/Tm_dG/max_dg_fraction
    if ss_dG >= 0:
        ss_score = 0
    if re.search(r'((\w)\2{5,})',seq):
        ss_score = 1e99
    tm_score = ((68-Tm)/68)**2    
    #score = (off_target_score + ss_score)/2
    score = 0.2*ss_score + 0.2*tm_score + 0.6*off_target_score + sum_target_score

    return score    

def score_q_consec(d,max_off_targets,maxid):
    # Extract variables from the dictionary
    ss_dG = float(d['ss_dG'])
    off_target_no = int(d['off_target_no'])
    Tm_dG = float(d['Tm_dG'])
    seq = str(d['sequence'])
    # Parameters
    max_dg_fraction = 0.5
    # Calculate the score
    if off_target_no > max_off_targets:
        off_target_score = 1e99
    else:
        off_target_score = off_target_no/max_off_targets
    ss_score = 1
    if(ss_dG/Tm_dG <= max_dg_fraction):
        ss_score = ss_dG/Tm_dG/max_dg_fraction
    if ss_dG >= 0:
        ss_score = 0
    query = r'((\w)\2{'+str(maxid-1)+',})'
    if re.search(query,seq):
        ss_score = 1e99

    score = (off_target_score + ss_score)/2
    return score

def score_q_combined(d,max_consec,maxid,targetTemp):         
   
    # Extract variables from the dictionary
    ss_dG = float(d['ss_dG'])
    off_target_no = int(d['off_target_no'])    # actually the longest consecutive stretch
    off_target_sum = int(d['off_target_sum'])  # cumulated mismatch counts
    Tm_dG = float(d['Tm_dG'])
    Tm = float(d['Tm'])
    seq = str(d['sequence'])

    # Parameters
    max_dg_fraction = 0.5

    # Calculate the oligo cost
    # 1. Secondary structures
    ss_cost = 1e10
    if ss_dG >= 0:
        ss_cost = 0
    elif(ss_dG/Tm_dG <= max_dg_fraction):
        ss_cost = ss_dG/Tm_dG/max_dg_fraction       # between 0-1 

    # 2. Melting temperature
    # Build into escafish in the future to minimize the Tm spread at the probe level?
    tm_cost = ((targetTemp-Tm)/5)**2  # +-5 deg should give at 0-1 range

    # 3. Off-target homology. Consecutive stretches
    if off_target_no > max_consec:
        consec_cost = 1e10
    else:
        consec_cost = off_target_no/max_consec     # between 0-1

    # 4. Off-target homology. Cumulated Hamming distances
    # The hard filter on consecutive stretches restrict the range of values to be expected here.
    sum_offtarget_cost = (10/(1+off_target_sum))**(1/2)     # between 0-1 (min value for a central mismatch in a 30mer: 9)

    # 5. Homopolymers
    query = r'((\w)\2{'+str(maxid-1)+',})'      
    if re.search(query,seq):
        ss_cost = 1e10

    # Regroup into single oligo cost
    score = 0.2*ss_cost + 0.2*tm_cost + 0.6*consec_cost + 0.2*sum_offtarget_cost

    return score    

def score_q_combined_bl(d,max_consec,maxid,targetTemp,hamdist):         
   
    # Extract variables from the dictionary
    ss_dG = float(d['ss_dG'])
    off_target_no = int(d['off_target_no'])    # actually the longest consecutive stretch
    off_target_sum = int(d['off_target_sum'])  # cumulated mismatch counts
    Tm_dG = float(d['Tm_dG'])
    Tm = float(d['Tm'])
    seq = str(d['sequence'])
    distBL = int(d['bl_dist'])

    # Parameters
    max_dg_fraction = 0.5

    # Calculate the oligo cost
    # 1. Secondary structures
    ss_cost = 1e10
    if ss_dG >= 0:
        ss_cost = 0
    elif(ss_dG/Tm_dG <= max_dg_fraction):
        ss_cost = ss_dG/Tm_dG/max_dg_fraction       # between 0-1 

    # 2. Melting temperature
    # Build into escafish in the future to minimize the Tm spread at the probe level?
    tm_cost = ((targetTemp-Tm)/5)**2  # +-5 deg should give at 0-1 range

    # 3. Off-target homology. Consecutive stretches
    if off_target_no > max_consec:
        consec_cost = 1e10
    else:
        consec_cost = off_target_no/max_consec     # between 0-1

    # 4. Off-target homology. Cumulated Hamming distances
    # The hard filter on consecutive stretches restrict the range of values to be expected here.
    if distBL <= hamdist:
        sum_offtarget_cost = 1e10
    else:
        sum_offtarget_cost = (10/(1+off_target_sum))**(1/2)     # between 0-1 (min value for a central mismatch in a 30mer: 9)

    # 5. Homopolymers
    query = r'((\w)\2{'+str(maxid-1)+',})'      
    if re.search(query,seq):
        ss_cost = 1e10

    # Regroup into single oligo cost
    score = 0.2*ss_cost + 0.2*tm_cost + 0.6*consec_cost + 0.2*sum_offtarget_cost

    return score    


score_functions = {'gg' : score_gg,
                   'gg_nhush': score_gg_nhush,
                   'q': score_q_nhush,
                   'q_cc': score_q_consec,
                   'q_combined': score_q_combined,
                   'q_bl': score_q_combined_bl}

if __name__ == '__main__':
    header = sys.stdin.readline().strip()

    if(len(sys.argv) < 2):
        print(f'No score function selected')
        print(f'Available:')
        for key in score_functions:
            print(f'{key}')
        exit(-1)

    score_function = score_functions[sys.argv[1]]

    # unless specified, allow max 6 consecutive identical base pairs
    maxid = "6"
    if(len(sys.argv) == 3):
        maxmm = int(sys.argv[2])
    elif(len(sys.argv) == 5):
        maxmm = int(sys.argv[2])
        maxid = int(sys.argv[3])
        targetTemp = int(sys.argv[4])
    elif(len(sys.argv) == 6):
        maxmm = int(sys.argv[2])
        maxid = int(sys.argv[3])
        targetTemp = int(sys.argv[4])    
        hamdist = int(sys.argv[5]) 
    print(f"{header}\toligo_cost")

    for line in tqdm(sys.stdin,"Reading database"):
        line = line.strip()
        d = dict(zip(header.split('\t'), line.split('\t')))

        if (sys.argv[1] == 'q_cc'):
            score = score_function(d,maxmm,maxid)
        elif (sys.argv[1] == 'q_combined'):
            score = score_function(d,maxmm,maxid,targetTemp)
        elif (sys.argv[1] == 'q_bl'):
            score = score_function(d,maxmm,maxid,targetTemp,hamdist)
        else: 
            score = score_function(d)
        # Optimization will not work if score < 0
        if(score < 0):
            score = 0

        print(f"{line}\t{score}")
