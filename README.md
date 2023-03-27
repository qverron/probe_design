# Instructions for probe design

CAUTION: These instructions are up to date!

## Preparation

### General:

- Install [nHUSH](https://github.com/elgw/nHUSH)

- Install [oldHUSH](https://github.com/elgw/hush)

- Install [escafish](https://github.com/elgw/escafish)

- Install [OligoArrayAux](http://www.unafold.org/Dinamelt/software/oligoarrayaux.php)

- Install [ifpd2](https://github.com/ggirelli/ifpd2)

### DNA:

- Get the genomic coordinates of the regions of interest

- Get the reference genome

### RNA:

- Get the transcripts of interest

- Get the reference transcriptome

### Notes:

- For DNA probes, the reference genome will be used both to extract
  the sequences of interest and to test probe candidates for
  homology. If different genomes need to be used, follow RNA steps and
  provide the regions of interest directly.
  
- For combined DNA-RNA FISH, the probe sets should be designed with an 
  homology check against both genome and transcriptome.
  
- All the commands below assume you are starting from your pipeline
  installation folder:
  ``` shell
  cd probe_design
  ```

## Probe design pipeline:

1. Preparation
- The probe desin pipeline data is currently intended to be run on a 
  folder called `data/` contained within the pipeline installation folder.

- Upon starting the pipeline, the `data/` folder should only contain
  `data/rois/` and `data/ref/`. If more folders are included, consider 
  making a back-up or simply removing them.
  
- List your regions of interest and their coordinates in the input file:
  `data/rois/all_regions.tsv`

- (Alternative: add  DNA / RNA separately  in `fw_DNA/RNA_FISH.txt` then
  run  `Rscript prepare_input.r`, which will gather them with correct
  formatting in `all_regions.tsv`)

- Place your reference genome in the `data/ref/` folder. Make
  sure that the chromosome naming matches with the reference genome
  name provided in `all_regions.tsv`.

- Generate all required subfolders:

  ``` shell
  mkdir data/candidates
  mkdir data/melt
  mkdir data/secs
  mkdir data/db
  mkdir data/db_tsv
  mkdir data/logfiles
  ```

2. Retrieve your region sequences and extract all k-mers of correct length:

   ``` shell
   # (from Pipeline/)
   ./get_oligos.py DNA|RNA [optional: applyGCfilter 0|1]
   # Example:
   ./get_oligos.py DNA 1
   ```

   If indicating RNA, the module will assume that the transcript / region
   sequences are already present in the data/regions folder. Default: DNA


3. Test all k-mers for their homology to other regions in the genome,
   using nHUSH. Instead of running the entire k-mers (of length `L`) at
   once, can be sped up by testing shorter sublength oligos (of length
   l).  `-m` number of mismatches to test for (always use 1 when running
   sublength); `-t` number of threads, `-i` comb size


- Full length:

  ``` shell
  ./run_nHUSH.sh -d RNA -L 35 -m 5 -t 40 -i 14
  ```
- Sublength:
  ``` shell
  ./run_nHUSH.sh -d DNA -L 40 -l 21 -m 1 -t 40 -i 14
  ```
- In case nHUSH is interrupted before completion, run before continuing:
  ``` shell
  ./unfinished_HUSH.sh
  ```
  
  4. Recapitulate nHUSH results as a score 

Recommended:

``` shell
./reform_hush_combined.py DNA|RNA|-RNA length sublength until
```
(`until` denotes the same number as specified after `-m` when running nHUSH).

- OR: If using sublength, recapitulate the mismatches into maximum
  consecutive match found:

  ``` shell
  ./reform_hush_consec.py DNA 40 21
  ```

- OR: When using either full-length or sublength (deprecated):
  ``` shell
  ./reform_hush.py DNA 40 21
  ```
  (syntax: `./reform_hush.py DNA|RNA|-RNA length (sublength)`)
  

5. Calculate the melting temperature of k-mers and the free energy of
   secondary structure formation:

   ``` shell
   ./melt_secs_parallel.sh (optional DNA(ref) / RNA(rev. compl))
   ```

6. Create k-mer database, convert to TSV for querying and attribute
   score to each oligo (based on nHUSH score, GC content, melting
   temperature, homopolymer stretches, secondary structures).
   
   Recommended:
   
   ``` shell
	./build-db.sh q_combined 32 6 70
    ```
    (optional score function: `q/gg/gg_nhush/q_combined`, default: `q`.
    Recommended: `q_combined`.

    32: Maximum length of a consecutive match. Default: 24
    6: Maximum length of a consecutive homopolymer. Default: 6
    All oligos with a longer consecutive match or homopolymer are stricly excluded.
    70: Target melting temperature. Default: 72C)
  
   Alternative:

    ``` shell
	./build-db_cc.sh q_cc 32
    ```
    (`q_cc` can only be used with sublength nHUSH and after running `reform_hush_consec.py`.

    32: Length of the maximum possible consecutive match.
	All oligos with a longer consecutive match are stricly excluded.)


7. (new!) Query the database to get candidate probes:

    ``` shell
	./cycling_query.py -s RNA -L 30 -m 7 -c 50 -t 40            
	[optional: -start 20 -end 100 -step 5]	 
To sweep different oligo numbers, otherwise uses the oligo counts provided in `./rois/all_regions.tsv`
        [optional: -stepdown 10]
Number of oligos to decrease probe size with every iteration that does not find enough oligos. Default: 1
    ```
Cycling query which generate probe candidates, then checks the resulting oligos using HUSH, removes inacceptable oligos and generate probes again.
If enough oligos cannot be found, design probes with fewer oligos, decreasing with `stepdown` at each step.

8. Summarize the final probes:

	If using `q_combined` scoring function:

   ``` shell
   ./summarize-probes-final.py
   ```

	If using one of the deprecated scoring functions:

   ``` shell
   python summarize-probes.py
   ```
	Some visual elements can be obtained using the following notebooks (needs updating!):
    ```shell
    plot_probe_candidates.ipynb
    plot_oligos.ipynb
    ```


## Generate probes for ordering

- Select forward, reverse primers and color flaps.
- Add the forward and reverse primer sequences to the probe oligos
- The forward primer to order has the color flap + the forward sequence
- The reverse primer to order has the t7 promoter sequence + the
  rev. compl of the rev sequence in the oligo
- The complete oligos can be uploaded as an Excel file containing the
  oligo names (arbitrary but unique) and the sequences
  
  
## TO DO:
- Adapt the code for more flexibility in input/output folders.
- Add a visual report of the probes at the end of the pipeline.
- One-button process!
- Find a way to automatize selecting primer sequences.

