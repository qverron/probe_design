set -e

datapath="$PWD/data"
suffix="Reference"
length=40
sublength=21

#if [ ! -z "$sublength" ]
#then
#	for d in "$datapath"/candidates/*"$suffix".fa."$sublength"mers
#        do
#            nhush dump-mindist "$d" "$d".mindist.uint8 "$sublength"
#		done
#else
#   for d in "$datapath"/candidates/*"$suffix".fa
#		do
#			nhush dump-mindist "$d" "$d".mindist.uint8 "$length"		
#		done
#fi


# reshape FASTA files in place
cd "$datapath"
for f in candidates/*"$suffix".fa; do echo $f; sed -r 's/pos=([0-9A-Za-z_]+):([0-9]+)-([0-9]+)\|([0-9]+):([0-9]+)/ \1 \2 \3 \4 \5/' $f | awk ' /^>/ {print $1" pos="$2":"$3+$5-1"-"$3+$6-1;next}1' > candidates/$(basename $f).fix; done
rm -r candidates/*"$suffix".fa
rename 's/.fix$//' candidates/*
