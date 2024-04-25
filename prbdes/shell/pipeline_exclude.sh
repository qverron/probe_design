# run all the steps of the pipeline
# 0. if necessary, retrieve the reference genome
./get_ref_genome.sh           # /!\ the settings need to be adjusted manually!!

# 0. Create all the directories for output
mkdir data/candidates
mkdir data/melt
mkdir data/secs
mkdir data/db 
mkdir data/db_tsv
mkdir data/probe_candidates
mkdir data/selected_probes

mkdir HUSH

# 1. Reformat input
Rscript prepare_input.r

# 2. Retrieve sequences and extract k-mers
./get_oligos.py DNA/RNA [applyGCfilter False] [extfolder] # RNA assumes already existing transcript sequences. Default> DNA

# 3. The exclusion pipeline requires a mask of regions to exclude from HUSH analysis. 
# These can either be added manually in data/exclude or automatically generated to only include the ROI itself. 
# Minimum: (use as template if adding more masks)
./generate_exclude.py

# 4. Generate masked reference genomes
./exclude_region.py 

# 5. Generate black list of abundantly repeated oligos for each reference genome
./generate_blacklist.sh -L 40 -c 100
# L: oligo length; c: min abundance to be included in oligo black list

# 6. Run nHUSH and reconstitute into full-length oligos if using sublength hashing
./run_nHUSH_excl.sh -d DNA -L 40 -l 21 -m 3 -t 40 -i 14 # 40 threads for max perf.

# only in case HUSH didn't finish successfully
./unfinished_HUSH.sh

# 7. only if using sublength:
./reform_hush_combined.py DNA|RNA|-RNA length sublength until

# 8. Calculate secondary structures and melting temperature
./melt_secs_parallel.sh DNA|RNA

# 9. Create database and convert to TSV for querying. Attribute score to each oligo
./build-db.sh q_combined 32 6 70 #(optional score function: q/q_combined, default: q)
./build-db_cc.sh q_cc 32 6  #score function: q_cc
./build-db_BL.sh -f q_bl -m 32 -i 6 -L 40 -c 100 -d 8 -T 72 # f: score function
# d: max Hamming distance to blacklist that is excluded
# L: oligo length; c: min abundance to be included in oligo black list
# i: max identical consecutive base pairs, T: target temperature, m: max length of consecutive off-target match

# 10. Query to fetch and optimize probes
./cycling_query.py -s DNA -L 40 -m 8 -c 100 -t 40 -g 500 -stepdown 50 -greedy -excl
# [recommended: -greedy     for speed > quality. Can be removed when designing final oligos]        
# [optional: -start 20 -end 100 -step 5: to sweep different oligo numbers]
# -stepdown 10       how many oligos to decrease with every iteration 

# 11. Summarize the final probes
python summarize-probes-final.py

# TO DO: ADD PLOTS
plot_oligos.ipynb
plot_probe_candidates.ipynb


