# run all the steps of the pipeline
# if necessary, retrieve the reference genome
# ./get_ref_genome.sh           # /!\ the settings need to be adjusted manually

# 0. Create all the directories for output
mkdir data/candidates
mkdir data/melt
mkdir data/secs
mkdir data/db 
mkdir data/db_tsv

mkdir HUSH

# 1. Reformat input
Rscript prepare_input.r

# 2. Retrieve sequences and extract k-mers
./get_oligos.py

# 3. Run nHUSH and reconstitute into full-length oligos if using sublength hashing
./run_nHUSH.sh -d DNA -L 40 -l 21 -m 2 -t 40 -i 14 # 40 threads for max perf. Add -g for genome ref
./run_nHUSH.sh -d RNA -L 35 -m 5 -t 40 -i 40

# only in case HUSH didn't finish successfully
./unfinished_HUSH.sh

# only if using sublength:
./reform_hush.py DNA 40 21  # syntax: ./reform_hush.py DNA|RNA|-RNA length sublength

# 4. Calculate secondary structures and melting temperature
./melt_secs_parallel.sh

# 5. Create database and convert to TSV for querying. Attribute score to each oligo
./build-db.sh

# 6. Query to fetch potential probe
./probe-query.sh -p 0.0001 -o 48

# 7. Inspect probes
python summarize-probes.py
# notebook
plot_probes.ipynb

