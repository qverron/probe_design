#!/bin/bash

##########################################

Help()
{
   # Display Help
   echo "Compare oligo database to blacklist."
   echo "Processes all databases in data/db_tsv," 
   echo "and fetches the blacklist files in data/blacklist."
   echo "If blacklists were generated for masked genomes,"
   echo "the database is compared to the corresponding blacklist."
   echo ""
   echo "Syntax: ./apply_blacklist.sh -L length -c abundance -d distance"
   echo "Arguments:"
   echo "L     kmer length"
   echo "c     Min abundance of an oligo to be included in the blacklist"
   echo "d     Minimum Hamming distance to any oligo in the blacklist."
   echo ""
   echo "Options:"
   echo "h     Show help"
}

##########################################
# Variables
length=40
cutoff=100
hamdist=9

while getopts "d:L:c:h" flag; do
   case "${flag}" in
      d) hamdist=${OPTARG};;
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
BLfolder="$data"/blacklist/
   
for dbfile in "$data"/db_tsv/db*.tsv
   do 
      roi=`echo $(basename -- "$dbfile") | sed 's/.*.\(roi_[0-9]\+\).*/\1\.fa/'`
      roiBL="$BLfolder"genome_"$roi".abundant_L"$length"_T"$cutoff".fa
      genomeBL="$BLfolder"genome.fa.abundant_L"$length"_T"$cutoff".fa
      if [ -f "$roiBL" ]
      then
         escafish apply_blacklist --db "$dbfile" --bl "$roiBL" 
      else
         escafish apply_blacklist --db "$dbfile" --bl "$genomeBL"   
      fi
done      
echo "Done!"