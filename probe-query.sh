#!/bin/bash
#set -e

#read -a arr < demo.txt 
#$ echo ${arr[0]}


Help()
{
   # Display Help
   echo "Construct optimized probes from pre-formed oligo databases."
   echo 
   echo "Syntax: probe-query -p 0.0001 -o 50"
   echo "-p pair weight (default 0.0001)"
   echo "If not provided, scanning pair weights between 1e-1 and 1e-7."
   echo "-o oligos per probe (default 50)"
   echo "If not provided, reading oligo number from ROI list."
   echo "-h display help"
}

while getopts "p:o:s:h" flag; do
   case "${flag}" in
      p) pw=${OPTARG};;
      o) oligos=${OPTARG};;
      s) type=${OPTARG};;
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
input="$expfolder/db_tsv"
suffix='Reference'
if [ ! -z "$type" ]
then 
case "$type" in
	DNA) suffix="Reference";;
	RNA) suffix="RevCompl";;
esac
fi

if [ -z $pw ]
then
for pwi in 1E-1 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7
do

ts="$(date +%Y%m%d-%H%M%S)"
randomid=$( cat /dev/urandom | tr -cd '0-9' | head -c 6 )       # to avoid folder names colliding
output="$expfolder/probe_candidates/query_output_t_$ts""-$randomid""_p_$pwi"

mkdir "$expfolder/probe_candidates"
mkdir $output

while IFS=$'\t' read -ra table
do
    echo "Processing region ""${table[2]}"
    roi="${table[2]}"
    file="$input/db.roi_$roi.GC35to85_$suffix.tsv"
    if [ -z $oligos ]
    then
        roioligos="${table[6]}"
    else
        roioligos="$oligos"   
    fi 
    out="$output/probe_roi_$roi.$roioligos""oligos.tsv"
    echo "Constructing probe with $roioligos oligos for region $roi"
    escafish --db $file --noligos $roioligos --out $out --pw $pwi --opt

done < <(tail -n +2 "$expfolder"/rois/all_regions.tsv)
done

else

ts="$(date +%Y%m%d-%H%M%S)"
randomid=$( cat /dev/urandom | tr -cd '0-9' | head -c 6 )       # to avoid folder names colliding
output="$expfolder/probe_candidates/query_output_t_$ts""-$randomid""_p_$pw"

mkdir $output

while IFS=$'\t' read -ra table
do
    echo "Processing region ""${table[2]}"
    roi="${table[2]}"
    file="$input/db.roi_$roi.GC35to85_$suffix.tsv"
    if [ -z $oligos ]
    then
        roioligos="${table[6]}"
    else
        roioligos="$oligos"   
    fi 
    out="$output/probe_roi_$roi.$roioligos""oligos.tsv"
    echo "Constructing probe with $roioligos for region $roi"
    escafish --db $file --noligos $roioligos --out $out --pw $pw --opt

done < <(tail -n +2 "$expfolder"/rois/all_regions.tsv)
fi



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
