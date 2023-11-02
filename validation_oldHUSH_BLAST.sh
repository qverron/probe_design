#!/bin/bash

############################################################

Help()
{
   # Display Help
   echo "Run HUSH and BLAST on the selected probe candidates."
   echo "The probe candidates should be located in their " 
   echo "respective query folders, and the reference genome in"
   echo "'./data/ref/genome.fa'."
   echo 
   echo "Syntax: validation_oldHUSH_BLAST -f folder -d FISH type -s direction -L length -m mism -t threads (-h|-p)"
   ech
   echo "Arguments:"
   echo "f     Experiment folder" 
   echo "d     FISH type (DNA|RNA)"
   echo "s     Strand containing the sequence of interest (p|n)"
   echo "L     kmer length"
   echo "m     Max number of mismatches being investigated"
   echo "t     Number of threads used for computing"
   echo "e     Exclude mode"
   echo
   echo "Options:"
   echo "h     Display help"
   echo "p     Skip division in sublength kmers"
}

##########################################
# Variables
gen=false
exclude=false

while getopts "f:L:m:t:he" flag; do
   case "${flag}" in
      f) exppath=${OPTARG};;
      L) length=${OPTARG};;
      m) mismatch=${OPTARG};;
      t) threads=${OPTARG};;
      e) exclude=true;;
      h) # display Help
         Help
         exit;;
     \?) # Invalid option
         echo "Error: Invalid option, exiting."
         exit;;
   esac
done

if [ -z "$exppath" ]
then
	exppath=$PWD
fi	
echo "Experiment folder: $exppath"
   
echo "Oligo length: $length"

echo "Mismatches: $mismatch"
echo "Threads: $threads"

if $gen
then
	echo "Assemble a genome reference"
fi	


# define a unique HUSH folder
datapath="${exppath}/data"
HUSHpath="${exppath}/HUSH"
ts="oldHUSH_"$(date +%Y%m%d-%H%M%S)
mkdir "$HUSHpath"/"$ts"

for gen in "$datapath"/ref/genome*.fa
do
ln $gen "$HUSHpath"/"$ts"/$(basename -- $gen)
done

# transform probe candidates into FASTA files readable by HUSH
for p in "$datapath"/selected_probes/*.tsv;
#do sed -r -n -e 's/^ROI_([0-9]+)'$'\t''chr([0-9A-Za-z_]+)'$'\t''([0-9]+)'$'\t''([0-9]+)'$'\t''([A-Z]+).*$/ROI_\1:\2:\3-\4|\5/p' $p | tr '|' '\n' > data/selected_probes/query_$(basename $p .tsv).fa;
do sed -r -n -e 's/^ROI_([0-9]+)'$'\t''([0-9A-Za-z_]+)'$'\t''([0-9]+)'$'\t''([0-9]+)'$'\t''([A-Z]+).*$/ROI_\1:\2:\3-\4|\5/p' $p | tr '|' '\n' > data/selected_probes/query_$(basename $p .tsv).fa;
done
wait

# run HUSH on the FASTA files
echo 'Creating subfolders'
for pfa in "$datapath"/selected_probes/*.fa;
do    
   queryfolder="$datapath"/selected_probes/$(basename $pfa .fa)_split_mm_$mismatch
   mkdir $queryfolder
done
wait   
echo 'Creating threads'
for pfa in "$datapath"/selected_probes/*.fa;
do  
   queryfolder="$datapath"/selected_probes/$(basename $pfa .fa)_split_mm_$mismatch
   cd $queryfolder
   # Create one input file per thread
   split -n l/$threads ../$(basename $pfa)
   cd ..
done
wait
echo 'Running HUSH'
for pfa in "$datapath"/selected_probes/*.fa;
do  
   if $exclude
   then
      genomeroi=`echo $(basename -- "$pfa") | sed 's/.*.\(roi_[0-9]\+\).*/genome_\1\.fa/'`
   else
      genomeroi="genome.fa"
   fi   
   queryfolder="$datapath"/selected_probes/$(basename $pfa .fa)_split_mm_$mismatch
   hushp -l $length -t $threads -r "$datapath"/ref/$genomeroi -q $queryfolder -m $mismatch -f 0 -C --verbose 1
done
wait
echo 'Exporting results'
for pfa in "$datapath"/selected_probes/*.fa;
do  
   queryfolder="$datapath"/selected_probes/$(basename $pfa .fa)_split_mm_$mismatch
   # Merge the ouputs to a single file
   cat $queryfolder/*.out > "$datapath"/selected_probes/$(basename $pfa)_"$mismatch"_mm.out
done
wait
echo 'Done!'


