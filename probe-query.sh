#!/bin/bash

#read -a arr < demo.txt 
#$ echo ${arr[0]}


Help()
{
   # Display Help
   echo "Construct optimized probes from pre-formed oligo databases."
   echo 
   echo "Syntax: probe-query -p 0.01 -o 50"
   echo "-p pair weight (default 0.01)"
   echo "-o oligos per probe (default 50)"
   echo "-h display help"
}

pw=0.01
oligos=50

while getopts "p:o:h" flag; do
   case "${flag}" in
      p) pw=${OPTARG};;
      o) oligos=${OPTARG};;
      h) # display Help
         Help
         exit;;
     \?) # Invalid option
         echo "Error: Invalid option, exiting."
         echo "probe-query -h for help"
         exit;;
   esac
done

expfolder="$PWD/data"
ts="$(date +%Y%m%d-%H%M%S)"
input="$expfolder/db_tsv"
output="$expfolder/query_output_$ts"

mkdir $output


while IFS=$'\t' read -ra table
do
    echo "Processing region ""${table[2]}"
    roi="${table[2]}"
    file="$input/db.roi_$roi.GC35to85_Reference.tsv"
    out="$output/probe_roi_$roi.$oligos""oligos.tsv"
    echo "Constructing probe with $oligos for region $roi"
    escafish --db $file --noligos $oligos --out $out --pw $pw

done < <(tail -n +2 "$expfolder"/rois/all_regions.tsv)



# # enable job control
# set -m

# process () {
#     local roi=$1
#     local file="$input/db.roi_$roi.GC35to85_Reference.tsv"
#     local out="$output/probe_roi_$roi.$oligos""oligos.tsv"
#     echo "Constructing probe with $oligos for region $roi"
#     probelet --db $file --noligos $oligos --out $out --pw $pw --iotn
# }

# for f in $(seq 1 $totalROI); do process $f & done
# # Wait for all parallel jobs to finish
# wait
# echo "Done!"