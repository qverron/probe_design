#!/bin/bash

############################################################

Help()
{
   # Display Help
   echo "Generate blacklist of abundantly repeated oligos in the genome."
   echo "Processes all ref genome files in data/ref," 
   echo "both the original and masked files."
   echo "'./data/ref/genome.fa'."
   echo 
   echo "Syntax: ./generate_blacklist.sh -L length -c abundance"
   echo "Arguments:"
   echo "L     kmer length"
   echo "c     Min abundance of an oligo to be included in the blacklist"
   echo ""
   echo "Options:"
   echo "h     Show help"
}

##########################################
# Variables
length=40
cutoff=100

while getopts "L:c:h" flag; do
   case "${flag}" in
      L) length=${OPTARG};;
      c) cutoff=${OPTARG};;
      h) # display Help
         Help
         exit;;
     \?) # Invalid option
         echo "Error: Invalid option, exiting."
         exit;;
   esac
done

data=$PWD'/data'

mkdir "$data"/blacklist
   
for genomefile in "$data"/ref/genome*.fa
   do 
      nhush find-abundant --file "$genomefile" --length "$length" --threshold "$cutoff" --out ./data/blacklist/$(basename -- "$genomefile").abundant_L"$length"_T"$cutoff".fa
done      
echo "Done!"