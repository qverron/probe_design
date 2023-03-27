#!/bin/bash

data=$PWD'/data'
cd $data
suffix='Reference'
if [ ! -z "$1" ]
then 
case "$1" in
	DNA) suffix="Reference";;
	RNA) suffix="RevCompl";;
esac
fi
mkdir melt
mkdir secs

# process all the ROI in parallel

# enable job control
set -m

process () {
    local file=$1
    cd $data
    echo "Calculating melting temperatures"
    melt_duplex -C -t DNA:DNA -o 0.05e-6 -n 1.04 -f 50 $file > "$data"/melt/$(basename $file) 2>/dev/null 
    cd "$data"/secs
    echo "Calculating secondary structures"
    hybrid-ss-min -n DNA -N 1.04 $file >/dev/null 2>&1
}

for f in "$data"/candidates/*"$suffix".fa; do process "$f" & done
# Wait for all parallel jobs to finish
wait
echo "Done!"

rename 's/.fa/.tsv/' "$data"/melt/*
