# run all the steps of the pipeline
# if necessary, retrieve the reference genome
# ./get_ref_genome      # /!\ the settings need to be adjusted manually

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
prb prepare_input

# 2. Generate black list of abundantly repeated oligos for each reference genome
prb generate_blacklist -L 40 -c 100
# only needs to be run once for a given reference genome

# 3. Retrieve sequences and extract k-mers
prb get_oligos DNA/RNA [applyGCfilter 0/1] [extfolder] # RNA assumes already existing transcript sequences. Default> DNA

# 4. Run nHUSH and reconstitute into full-length oligos if using sublength hashing
prb run_nHUSH -d DNA -L 40 -l 21 -m 3 -t 40 -i 14 # 40 threads for max perf. Add -g for genome ref
prb run_nHUSH -d RNA -L 35 -m 5 -t 40 -i 14

# only in case HUSH didn't finish successfully
prb unfinished_HUSH 

# 5. Recapitulate HUSH results
prb reform_hush_combined DNA|RNA|-RNA length sublength until

# 6. Calculate secondary structures and melting temperature
prb melt_secs_parallel DNA|RNA

# 7. Create database and convert to TSV for querying. Attribute score to each oligo
prb build-db q_combined 32 6 70 #(optional score function: q/q_combined, default: q)
prb build-db_BL -f q_bl -m 32 -i 6 -L 40 -c 100 -d 8 -T 72 #score function: q_cc
# 32: length of max consecutive perfect match allowed; 6: max number of consecutive identical base pairs, 70: target temperature

# 8. Query to fetch and optimize probes
prb cycling_query -s DNA -L 40 -m 8 -c 100 -t 40 -g 2000 -greedy -stepdown 10
# [recommended: -greedy   for speed > quality. Can be removed when designing final oligos]    
# [optional: -start 20 -end 100 -step 5: to sweep different oligo numbers]
# -stepdown 10    how many oligos to decrease with every iteration 
# -opt        Get the best probes, but takes a lot longer. Can be done after a quick iteration to remove poor oligos


# 7. Summarize the final probes
prb summarize-probes-final 

# TO DO: ADD PLOTS
prb plot_oligos
prb plot_probe_candidates


