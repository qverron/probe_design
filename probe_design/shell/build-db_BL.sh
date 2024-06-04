#!/bin/bash

##########################################

Help()
{
   # Display Help
   echo "Generate oligo database for each ROI."
   echo "Compare the oligos to the blacklist of oligos"
   echo "with the specified parameters."
   echo "Attribute score to each oligo using the"
   echo "specified cost function."
   echo ""
   echo "Syntax: ./build-db_BL.sh -L length -c abundance -d distance -f function"
   echo "-m max consec -i max homopolymer -T temperature"
   echo "Arguments:"
   echo "L     kmer length"
   echo "c     Min abundance of an oligo to be included in the blacklist"
   echo "d     Minimum Hamming distance to any oligo in the blacklist (default: 8)"
   echo "f     Escafish cost function. Default & recommended: q_bl"
   echo "m     Longest consecutive match allowed (default: 24 nt)"
   echo "i     Longest homopolymer allower (default: 6 nt)"
   echo "T     Target melting temperature (default: 72C)"
   echo ""
   echo "Options:"
   echo "h     Show help"
}

##########################################
# Variables
length=40
cutoff=100
hamdist=8
scoref="q_bl"
maxconsec=24
maxid=6
targetTemp=72

confirm=false

while getopts "d:L:c:f:m:i:T:hy" flag; do
   case "${flag}" in
      d) hamdist=${OPTARG};;
      L) length=${OPTARG};;
      c) cutoff=${OPTARG};;
      f) scoref=${OPTARG};;
      m) maxconsec=${OPTARG};;
      i) maxid=${OPTARG};;
      T) targetTemp=${OPTARG};;
      h) # display Help
         Help
         exit;;
      y) confirm=true;;
     \?) # Invalid option
         echo "Error: Invalid option, exiting."
         exit;;
   esac
done

data=$PWD'/data'
cd "$data"
if [ -d "db" ]
then 
   echo 'The directories need to be cleared to continue.'
   if [ "$confirm" = false ] ; then
      while true; do
            read -p "Do you wish to continue? (y/n)" yn
            case $yn in
                [Yy]* ) break;;
                [Nn]* ) exit;;
                * ) echo "Please answer yes (y) or no (n).";;
            esac
        done
    fi

   rm -r db
   rm -r db_tsv
   rm -r db_temp1
   rm -r db_temp2
fi

mkdir db
BLfolder="blacklist/"

rename 's/.fa$/.hush.out/' HUSH_candidates/*
echo "Building databases."
for f in HUSH_candidates/*; do f2=$(basename $f | sed 's/.hush.out$//'); f3=$(echo "$f2"|sed 's/^sequences_//'); ifpd2 db make -O HUSH_candidates/"$f2".hush.out -T melt/"$f2".tsv -S secs/"$f2".fa.ct db/db."$f3"; done

mkdir db_temp1
echo "Converting to TSV."
for d in db/*; do ifpd2 db dump $d > db_temp1/$(basename $d).tsv; done

mkdir db_temp2
mkdir db_tsv
for d in db_temp1/*; do echo $d; sed -r 's/'$'\t''111([0-9]+)987([0-9]+)'$'\t''/'$'\t''\1'$'\t''\2'$'\t''/' $d | sed -r $'s/off_target_no\t/off_target_no\toff_target_sum\t/' > db_temp2/$(basename $d); done

echo "Comparing oligo database to blacklist"
for dbfile in db_temp2/db*.tsv
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

echo "Attributing oligo score."
echo "Using the following score function: $scoref"
#for d in db_temp1/*; do cat $d | ../escafish_score.py "$scoref" "$2" "$3" > db_tsv/$(basename $d); done
for d in db_temp2/*_filtered.fa; do cat $d | prb escafish_score "$scoref" "$maxconsec" "$maxid" "$targetTemp" "$hamdist" > db_tsv/$(basename $d ".bl_filtered.fa"); done
rm -r db_temp1
rm -r db_temp2

cd ..
