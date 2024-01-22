# Instructions for probe design

> [!CAUTION]
> You may want to add python and pip as aliases for python3 and pip3

On your terminal;

```shell
echo "alias python='python3'" >> ~/.bashrc
echo "alias pip='pip3'" >> ~/.bashrc
source ~/.bashrc
```

## Installation

- Install **probe_design**  (also installs **ifpd2q**)

```shell
pip install probe_design
```

Add (pip) installation directory for the up-coming use cases;

On your terminal;

```shell
echo "export PRB=$(pip show probe_design | awk '/Location:/{print $2"/probe_design"}')" >> ~/.bashrc
source ~/.bashrc
```

- Install **oligo-melting**

On your terminal;

```shell
git clone http://github.com/ggirelli/oligo-melting ~/oligo-melting
cd ~/oligo-melting
pip install .
```

> [!NOTE]
> nHUSH, HUSH and escafish are private repositories

- Install the **dev branch** of [nHUSH](https://github.com/elgw/nHUSH/tree/dev)

- Install [HUSH](https://github.com/elgw/hush)

- Install [escafish](https://github.com/elgw/escafish)

- Install [OligoArrayAux](http://www.unafold.org/Dinamelt/software/oligoarrayaux.php)


## Preparation

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
  
- All the commands below assume you are starting from your project directory

```shell
mkdir <project_name>
cd <project_name>
```

# Probe design pipeline:

## Alternative 1: Normally repetitive regions.

1. Preparation

- The probe desin pipeline data is currently intended to be run on a 
  folder called `data/` contained within the pipeline installation folder.

- Upon starting the pipeline, the `data/` folder should only contain
  `data/rois/` and `data/ref/` (and possibly `data/blacklist/`, see 6.). If more folders are included, consider making a back-up or simply removing them.
  
- List your regions of interest and their coordinates in the input file:
  `data/rois/all_regions.tsv`

- Place your reference genome in the `data/ref/` folder. Make
  sure that the chromosome naming matches with the reference genome
  name provided in `all_regions.tsv`.
  
- The reference folder can alternatively be gathered using `get_ref_genome.sh`.
  In that case, adjust the script manually with the correct Ensembl 
  address for your genome of interest.

2. Generate all required subfolders inside your project directory:

```shell
mkdir data/candidates
mkdir data/melt
mkdir data/secs
mkdir data/db
mkdir data/db_tsv
mkdir data/logfiles
mkdir HUSH
```  

3. Retrieve your region sequences and extract all k-mers of correct length:

```shell
# ($PRB is the path to pip installation of probe_design)
# src is for the python source codes
python $PRB/src/get_oligos.py DNA|RNA [optional: applyGCfilter 0|1]
# Example:
python $PRB/src/get_oligos.py DNA 1
```

> [!NOTE]
> If indicating `RNA`, the module will assume that the transcript / region
> sequences are already present in the `data/regions` folder. Default: `DNA.


4. Test all k-mers for their homology to other regions in the genome,
   using nHUSH. Instead of running the entire k-mers (of length `L`) at
   once, can be sped up by testing shorter sublength oligos (of length
   l).  `-m` number of mismatches to test for (always use 1 when running
   sublength); `-t` number of threads, `-i` comb size

- Full length:

``` shell
bash $PRB/shell/run_nHUSH.sh -d RNA -L 35 -m 5 -t 40 -i 14
```

- Sublength:

``` shell
bash $PRB/shell/run_nHUSH.sh -d DNA -L 40 -l 21 -m 3 -t 40 -i 14
```

> [!TIP]
> ADD -g if this is the first time running with a new reference genome!  
  
- In case nHUSH is interrupted before completion, run before continuing:

``` shell
bash $PRB/shell/unfinished_HUSH.sh
```
  
5. Recapitulate nHUSH results as a score 

``` shell
python $PRB/src/reform_hush_combined.py DNA|RNA|-RNA length sublength until
```

(`until` denotes the same number as specified after `-m` when running nHUSH). 

6. Calculate the melting temperature of k-mers and the free energy of
   secondary structure formation:

``` shell
bash $PRB/shell/melt_secs_parallel.sh (optional DNA(ref) / RNA(rev. compl))   
```

7. Generate a black list of abundantly repeated oligos in the reference genome.

``` shell
bash $PRB/shell/generate_blacklist.sh -L 40 -c 100
```
> [!NOTE]
> This only needs to be run once per reference genome if not using any 
> exclusion regions! Just save the blacklist folder between runs.

> L: oligo length <br>
> c: min abundance to be included in oligo black list

8. Create k-mer database, convert to TSV for querying and attribute
   score to each oligo (based on nHUSH score, GC content, melting
   temperature, homopolymer stretches, secondary structures).

``` shell
bash $PRB/shell/build-db_BL.sh -f q_bl -m 32 -i 6 -L 40 -c 100 -d 8 -T 72
```

> m: Maximum length of a consecutive match. Default: 24 <br>
> i: Maximum length of a consecutive homopolymer. Default: 6 <br>
> All oligos with a longer consecutive match or homopolymer are stricly excluded. <br>
> L: oligo length <br>
> c: min number of occurrences for an oligo to be counted in black list <br>
> (should match settings used in 6.) <br>
> d: min Hamming distance to an oligo in the blacklist for exclusion  <br>
> T: Target melting temperature. Default: 72C



9. Query the database to get candidate probes:

``` shell
python $PRB/src/cycling_query.py -s DNA -L 40 -m 8 -c 100 -t 40 -greedy
```

**[optional: -greedy. Speed > quality]
[optional: -start 20 -end 100 -step 5]**

To sweep different oligo numbers, otherwise uses the oligo counts provided in `./rois/all_regions.tsv`
        [optional: -stepdown 10]
Number of oligos to decrease probe size with every iteration that does not find enough oligos. Default: 1

Cycling query which generate probe candidates, then checks the resulting oligos using HUSH, removes inacceptable oligos and generate probes again.
If enough oligos cannot be found, design probes with fewer oligos, decreasing with `stepdown` at each step.

10. Summarize the final probes:

```shell
python $PRB/src/summarize-probes-final.py
```

Some visual elements can be obtained using the following notebooks (needs updating!):

``` shell
python $PRB/notebooks/plot_probe_candidates.ipynb
python $PRB/notebooks/plot_oligos.ipynb
```


## Alternative 2: Repetitive or repeated regions.

In this alternative, the region (along with any user-indicated repeats)
is masked out from the reference genome used by nHUSH. This way, repeated
oligos that are specific for the ROI can be included in the final probe.

### Warning: This approach occupies a lot more hard drive space!

1. Preparation
- Besides `data/rois/` and `data/ref/`, the pipeline requires an additional
  `data/exclude/` folder containing BED files with the coordinates of sections
  to mask out when running HUSH for each ROI. 
  
2. (UNLESS manually providing exclusion regions)
Exclude regions of interest from HUSH scan.

``` shell
python $PRB/src/generate_exclude.py
```
- The same sheet template can be used to manually add further regions to exclude.

2. Generate all required subfolders:

``` shell
mkdir data/candidates
mkdir data/melt
mkdir data/secs
mkdir data/db
mkdir data/db_tsv
mkdir data/logfiles
```

3. Retrieve your region sequences and extract all k-mers of correct length:

``` shell
# (from Pipeline/)
python $PRB/src/get_oligos.py DNA|RNA [optional: applyGCfilter 0|1]
# Example:
python $PRB/src/get_oligos.py DNA
```

   If indicating `RNA`, the module will assume that the transcript / region
   sequences are already present in the `data/regions` folder. Default: `DNA.
   
4. Apply the region exclusion mask on the reference genome.

``` shell
python $PRB/src/exclude_region.py 
```   
   
5. Generate a black list of abundantly repeated oligos in the reference genome.
	
``` shell
bash $PRB/shell/generate_blacklist.sh -L 40 -c 100
```
Needs to be re-run everytime when using exclusion masks.
L: oligo length; c: min abundance to be included in oligo black list   


6. Test all k-mers for their homology to other regions in the genome,
   using nHUSH. Instead of running the entire k-mers (of length `L`) at
   once, can be sped up by testing shorter sublength oligos (of length
   l).  `-m` number of mismatches to test for (minimum 1 for sublength;
   more gives better information but takes longer time);
   `-t` number of threads, `-i` comb size

Sublength:

``` shell
bash $PRB/shell/run_nHUSH_excl.sh -d DNA -L 40 -l 21 -m 3 -t 40 -i 14
```
  
Note the `_excl` specific to the exclusion mode.  
  
In case nHUSH is interrupted before completion, run before continuing:
``` shell
bash $PRB/shell/unfinished_HUSH.sh
```
  
7. Recapitulate nHUSH results as a score 

``` shell
# Format:
python $PRB/src/reform_hush_combined.py DNA|RNA|-RNA length sublength until
# Example:
python $PRB/src/reform_hush_combined.py DNA 40 21 3
```

(`until` denotes the same number as specified after `-m` when running nHUSH).

8. Calculate the melting temperature of k-mers and the free energy of
   secondary structure formation:

``` shell
bash $PRB/shell/melt_secs_parallel.sh (optional DNA(ref) / RNA(rev. compl))
```

9. Create k-mer database, convert to TSV for querying and attribute
   score to each oligo (based on nHUSH score, GC content, melting
   temperature, homopolymer stretches, secondary structures).
   
   Recommended:
   
``` shell
bash $PRB/shell/build-db_BL.sh -f q_bl -m 32 -i 6 -L 40 -c 100 -d 8 -T 72
```
    
> f: score function  <br>
> d: max Hamming distance to blacklist that is excluded <br>
> L: oligo length <br>
> c: min abundance to be included in oligo blacklist <br>
> i: max identical consecutive base pairs,  <br>
> T: target temperature <br>
> m: max length of consecutive off-target match <br>
  
10. Query the database to get candidate probes:

``` shell
python $PRB/src/cycling_query.py -s DNA -L 40 -m 8 -c 100 -t 40 -g 500 -stepdown 50 -greedy -excl
```

**[optional: -greedy. Speed > quality]
[optional: -start 20 -end 100 -step 5]**

To sweep different oligo numbers, otherwise uses the oligo counts provided in `./rois/all_regions.tsv`
        [optional: -stepdown 10]
Number of oligos to decrease probe size with every iteration that does not find enough oligos. Default: 1

Cycling query which generate probe candidates, then checks the resulting oligos using HUSH, removes inacceptable oligos and generate probes again.
If enough oligos cannot be found, design probes with fewer oligos, decreasing with `stepdown` at each step.


11. Summarize the final probes:

``` shell
python $PRB/src/summarize-probes-final.py
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
