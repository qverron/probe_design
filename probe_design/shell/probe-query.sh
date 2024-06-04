#!/bin/bash
set -m      # enable job control


Help()
{
   # Display Help
   echo "Construct optimized probes from pre-formed oligo databases."
   echo 
   echo "Syntax: probe-query -L 40 -p 0.0001 -o 50"
   echo "-L oligo length"
   echo "-p pair weight (default 0.0001)"
   echo "If not provided, scanning pair weights between 1e-1 and 1e-7."
   echo "-o oligos per probe (default 50)"
   echo "If not provided, reading oligo number from ROI list."
   echo "-e specific ROI to query"
   echo "If not provided, querying all ROIs listed in all_regions.tsv"
   echo "-g Greedy probe query, speed > quality"
   echo "-h display help"
}

greedy=''

while getopts "L:p:o:s:e:l:hg" flag; do
   case "${flag}" in
      L) length=${OPTARG};;
      p) pw=${OPTARG};;
      o) oligos=${OPTARG};;
      s) type=${OPTARG};;
      e) probe=${OPTARG};;
      g) greedy='--greedy';;
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


ts="$(date +%Y%m%d-%H%M%S)"
randomid=$( cat /dev/urandom | tr -cd '0-9' | head -c 6 )       # to avoid folder names colliding
#output="$expfolder/probe_candidates/query_output_t_$ts""-$randomid""_p_$pwi"
if [ -z $oligos ]
then
    output="$expfolder/probe_candidates/query_output_t_$ts"   
else
    output="$expfolder/probe_candidates/query_output_o_$oligos""_t_$ts"
fi

mkdir -p "$expfolder/probe_candidates"
mkdir -p $output


process () {
    local pwl=$1
    while IFS=$'\t' read -ra table
    do
        roi="${table[2]}"
        if [ ! -z "$probe" ] && [ "$roi" != $probe ]        # the user has specified a probe, skip other regions
        then
            continue
        else    
            file="$input/db.roi_$roi.GC35to85_$suffix.tsv.filt"
            if [ -z $oligos ]
            then
                roioligos="${table[6]}"
            else
                roioligos="$oligos"   
            fi 
            out="$output/probe_roi_$roi.$roioligos""oligos.pw$pwl.tsv"
            echo "Constructing probe with $roioligos oligos for region $roi, pair weight: $pwl."
            #escafish --db $file --len $length --noligos $roioligos --out $out --pw $pwl $greedy
            escafish --db $file --noligos $roioligos --out $out --pw $pwl $greedy

        fi

    done < <(tail -n +2 "$expfolder"/rois/all_regions.tsv)
}


if [ -z $pw ]           # no pair weight specify, sweep over full range
then
for pwi in 1E-1 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7; do process "$pwi" & done
else        # single pair weight
process "$pw"
fi
# Wait for all parallel jobs to finish
wait
echo "Done!"
