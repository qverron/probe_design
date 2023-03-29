# run all the steps of the pipeline
# if necessary, retrieve the reference genome
# ./get_ref_genome.sh           # /!\ the settings need to be adjusted manually

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
./get_oligos.py DNA/RNA [applyGCfilter 0/1] [extfolder] # RNA assumes already existing transcript sequences. Default> DNA

# 3. Run nHUSH and reconstitute into full-length oligos if using sublength hashing
./run_nHUSH.sh -d DNA -L 40 -l 21 -m 1 -t 40 -i 14 # 40 threads for max perf. Add -g for genome ref
./run_nHUSH.sh -d RNA -L 35 -m 5 -t 40 -i 14

# only in case HUSH didn't finish successfully
./unfinished_HUSH.sh

# only if using sublength:
./reform_hush_consec.py DNA 40 21  # syntax: ./reform_hush.py DNA|RNA|-RNA length (sublength)
./reform_hush.py DNA|RNA|-RNA length (sublength)
./reform_hush_combined.py DNA|RNA|-RNA length sublength until

# 4. Calculate secondary structures and melting temperature
./melt_secs_parallel.sh DNA|RNA

# 5. Create database and convert to TSV for querying. Attribute score to each oligo
./build-db.sh q_combined 32 6 70 #(optional score function: q/q_combined, default: q)
./build-db_cc.sh q_cc 32 6  #score function: q_cc
# 32: length of max consecutive perfect match allowed; 6: max number of consecutive identical base pairs, 70: target temperature

# 6. Query to fetch and optimize probes
./cycling_query.py -s RNA -L 30 -m 7 -c 50 -t 40            
# [optional: -start 20 -end 100 -step 5: to sweep different oligo numbers]
# -stepdown 10       how many oligos to decrease with every iteration 
# -opt               Get the best probes, but takes a lot longer. Can be done after a quick iteration to remove poor oligos


# 7. Summarize the final probes
python summarize-probes-final.py

# TO DO: ADD PLOTS
plot_oligos.ipynb
plot_probe_candidates.ipynb


