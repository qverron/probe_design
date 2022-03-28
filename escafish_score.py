#!/usr/bin/python

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
        sum_target_score = 5/(1+off_target_sum)**(1/2)      
    ss_score = 1
    if(ss_dG/Tm_dG <= max_dg_fraction):
        ss_score = ss_dG/Tm_dG/max_dg_fraction
    if ss_dG >= 0:
        ss_score = 0
    if re.match(r'((\w)\2{5,})',seq):
        ss_score = 1e99
    tm_score = ((68-Tm)/68)**2    
    #score = (off_target_score + ss_score)/2
    score = 0.2*ss_score + 0.2*tm_score + 0.6*off_target_score + sum_target_score

    return score    

score_functions = {'gg' : score_gg,
                   'gg_nhush': score_gg_nhush,
                   'q': score_q_nhush}

if __name__ == '__main__':
    header = sys.stdin.readline().strip()

    if(len(sys.argv) < 2):
        print(f'No score function selected')
        print(f'Available:')
        for key in score_functions:
            print(f'{key}')
        exit(-1)

    score_function = score_functions[sys.argv[1]]
    print(f"{header}\toligo_cost")

    for line in tqdm(sys.stdin,"Reading database"):
        line = line.strip()
        d = dict(zip(header.split('\t'), line.split('\t')))

        score = score_function(d)
        # Optimization will not work if score < 0
        if(score < 0):
            score = 0

        print(f"{line}\t{score}")
