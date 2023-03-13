# Probe design

Autocorrected to menu
Conversation opened. 2 messages. 1 message unread.

Skip to content
Using Science for Life Laboratory Mail with screen readers
6 of 764
(no subject)
Inbox
Quentin Verron
	
	AttachmentsFri, Mar 10, 12:59 PM (3 days ago)
-- Quentin Verron Post-doctoral Researcher BiCro lab - Karolinska Institute / SciLifeLab Tomtebodavägen 23A - 17165 Solna https://bienkocrosettolabs.org/
Erik Wernersson
	
AttachmentsFri, Mar 10, 1:24 PM (3 days ago)
	
to me
Typ såhär.

Kanske inte ska ha med lösenordet om repot är publikt ...

E

 One attachment  •  Scanned by Gmail
	

# Instructions for probe design

CAUTION: These instructions are not up to date.

## Preparation

### DNA:
- Get the genomic coordinates of the regions of interest

- Get the reference genome. Commonly used reference genomes are
  already present in Probe_design/ref/ (see below). For other genomes,
  Pipeline/get_ref_genome.sh can be used but requires manual set-up.

### RNA:

- Get the transcripts of interest
- Get the reference genome (idem as for DNA)

### General:

- Set up SciLifeLab VPN, instructions here. You can check that the
  installation worked by successfully accessing the SciLifeLab
  intranet.

- Get familiar with basic Terminal commands, including how to navigate
  between folders, how to move or edit text files, etc.

### Notes:

- For DNA probes, the reference genome will be used both to extract
  the sequences of interest and to test probe candidates for
  homology. If different genomes need to be used, follow RNA steps and
  provide the regions of interest directly.



## Probe design pipeline:

1. Computer access
- Remotely access the computer used for probe design through SSH:
  ``` shell
  ssh probedesign@ggirellipc.scilifelab.se
  ```
  Password: probe2023
- The probe design folder can be found under /mnt/data/Probe_design/
- All pipeline functions are in Probe_design/Pipeline/
  and all input-output in Probe_design/Pipeline/data/
- Note: A step-by-step summary of the modules in the probe design
  pipeline is found in `Probe_design/Pipeline/pipeline.sh`


2. Preparation
- Check that no one else is running probe design!

- If `Probe_design/Pipeline/data/` contains data, rename the folder to
  `data_backup_initials_date` so the other person gets a chance to
  retrieve their data.

- Upon starting the pipeline, the data/ folder should only contain
  `data/rois/` and `data/ref/`.
- List your regions of interest and their coordinates in the input file:
  `Pipeline/data/rois/all_regions.tsv`

- (Alternative: add  DNA / RNA separately  in `fw_DNA/RNA_FISH.txt` then
  run  `Rscript prepare_input.r`, which will gather them with correct
  formatting in `all_regions.tsv`)

- Place your reference genome in the `Pipeline/data/ref/` folder. Make
  sure that the chromosome naming matches with the reference genome
  name provided in `all_regions.tsv`.

- Generate all required subfolders:

  ``` shell
  mkdir data/candidates
  mkdir data/melt
  mkdir data/secs
  mkdir data/db
  mkdir data/db_tsv
  ```

3. Retrieve your region sequences and extract all k-mers of correct length:

   ``` shell
   # (from Pipeline/)
   ./get_oligos.py (optional DNA/RNA)
   ```

   If indicating RNA, the module will assume that the transcript / region
   sequences are already present in the data/regions folder. Default: DNA


4. Test all k-mers for their homology to other regions in the genome,
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
- Recapitulate nHUSH results as a score
  ``` shell
  ./reform_hush.py DNA 40 21
  ```
  (syntax: `./reform_hush.py DNA|RNA|-RNA length (sublength)`)
- OR: If using sublength, recapitulate the mismatches into maximum
  consecutive match found:

  ``` shell
  ./reform_hush_consec.py DNA 40 21
  ```

5. Calculate the melting temperature of k-mers and the free energy of
   secondary structure formation:

   ``` shell
   ./melt_secs_parallel.sh (optional DNA(ref) / RNA(rev. compl))
   ```

6. Create k-mer database, convert to TSV for querying and attribute
   score to each oligo (based on nHUSH score, GC content, melting
   temperature, homopolymer stretches, secondary structures).

    ``` shell
	./build-db.sh q_cc 32
    ```
    (optional score function: `q/gg/gg_nhush/q_cc`, default: `q`.
    `q_cc` can only be used with sublength nHUSH and after running `reform_hush_consec.py`.

    32: Length of the maximum potential match allowed.)


7. Query the database to get candidate probes:

    ``` shell
	./probe-query.sh -s DNA/RNA
    ```
	Optional: `-o 48` (number of oligos to include in the probe)
    `-p 0.0001` (pair weight, otherwise sweeps range 1e-2 - 1e-7)

8. Inspect the generated probes:
   ``` shell
   python summarize-probes.py
   ```
	Some visual elements can be obtained using the following notebook:
    ```shell
    plot_probes.ipynb
    ```

9. Select a probe among the candidates.

## Generate probes for ordering

- Select forward, reverse primers and color flaps.
- Add the forward and reverse primer sequences to the probe oligos
- The forward primer to order has the color flap + the forward sequence
- The reverse primer to order has the t7 promoter sequence + the
  rev. compl of the rev sequence in the oligo
- The complete oligos can be uploaded as an Excel file containing the
  oligo names (arbitrary but unique) and the sequences

README.md
Displaying README.md.
