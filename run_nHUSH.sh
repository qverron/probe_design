#!/bin/bash

############################################################

Help()
{
   # Display Help
   echo "Run nHUSH on the selected probe candidates."
   echo "The probe candidates should be located in" 
   echo "'./data/candidates/' and the reference genome in"
   echo "'./data/ref/genome.fa'."
   echo 
   echo "Syntax: run_nHUSH -f folder -d FISH type -s direction -L length (-l sublength) -m mism -t threads -i hash (-h|-c|-g|-p)"
   ech
   echo "Arguments:"
   echo "f     Experiment folder" 
   echo "d     FISH type (DNA|RNA)"
   echo "s     Strand containing the sequence of interest (p|n)"
   echo "L     kmer length"
   echo "l     Sublength used for faster nHUSH (optional)"
   echo "m     Max number of mismatches being investigated"
   echo "t     Number of threads used for computing"
   echo "i     Initial hash length"
   echo
   echo "Options:"
   echo "h     Display help"
   echo "g     Assemble a reference genome from separate files"
   echo "p     Skip division in sublength kmers"
}

##########################################
# Variables
skip=false
gen=false

while getopts "f:L:l:m:t:i:d:s:hgcp" flag; do
   case "${flag}" in
      f) exppath=${OPTARG};;
      d) fishtype=${OPTARG};;
      s) strand=${OPTARG};;   
      L) length=${OPTARG};;
      l) sublength=${OPTARG};;
      m) mismatch=${OPTARG};;
      t) threads=${OPTARG};;
      i) inhash=${OPTARG};;
      h) # display Help
         Help
         exit;;
      g) gen=true;;
      p) skip=true;; 
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

case "$fishtype" in
	DNA) echo "DNA FISH"
	     suffix="Reference";;
	RNA) echo "RNA FISH"
             case "$strand" in
	     	p) echo "Positive strand"; suffix="RevCompl";;
	        n) echo "Negative strand"; suffix="Reference";;
		*) echo "Missing strand information, exiting."; exit;;
	     esac;;
esac	     

echo "Oligo length: $length"

if [ ! -z "$sublength" ]
then
	echo "Sublength: $sublength"
fi

if $skip
then echo "Skipping division in sublength kmers"
fi

echo "Mismatches: $mismatch"
echo "Threads: $threads"
echo "Hash length: $inhash"

if $gen
then
	echo "Assemble a genome reference"
fi	


# Check that the input were correctly parsed
while true; do
    read -p "Do you wish to continue?" yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

# define a unique nHUSH folder
datapath="${exppath}/data"
HUSHpath="${exppath}/HUSH"
ts="HUSH_"$(date +%Y%m%d-%H%M%S)
mkdir "$HUSHpath"/"$ts"

# assemble a complete reference genome
if $gen
then
	cat "$datapath"/ref/*.fa > "$datapath"/ref/genome.fa
fi

ln "$datapath"/ref/genome.fa "$HUSHpath"/"$ts"/genome.fa

if [ ! -z "$sublength" ]
then
   if ! $skip
	then
		for d in "$datapath"/candidates/*"$suffix".fa; do nhush-fasplit -L "$length" -l "$sublength" --file "$d"; done
	fi
	for d in "$datapath"/candidates/*"$suffix".fa."$sublength"mers
      do
         cd "$HUSHpath"/"$ts"
			nhush --hash "$inhash" --length "$sublength" --until "$mismatch" --threads "$threads" --external "$d" --file genome.fa --sfp
			#nhush dump-mindist "$d" "$d".mindist.uint8 "$sublength"
		done
else
   for d in "$datapath"/candidates/*"$suffix".fa
		do
			cd "$HUSHpath"/"$ts"
			nhush --hash "$inhash" --length "$length" --until "$mismatch" --threads "$threads" --external "$d" --file genome.fa --sfp
			#nhush dump-mindist "$d" "$d".mindist.uint8 "$length"		
		done
fi

# reshape FASTA files in place
cd "$datapath"
for f in candidates/*"$suffix".fa; do echo $f; sed -r 's/pos=chr([0-9X]+):([0-9]+)-([0-9]+)\|([0-9]+):([0-9]+)/ \1 \2 \3 \4 \5/' $f | awk ' /^>/ {print $1" pos=chr"$2":"$3+$5-1"-"$3+$6-1;next}1' > candidates/$(basename $f).fix; done
rm -r candidates/*"$suffix".fa
rename 's/.fix$//' candidates/*
