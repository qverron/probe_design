#!/bin/bash

cd data
if [ -d "db" ]
then 
    echo 'The directories need to be cleared to continue.'
    while true; do
        read -p "Do you wish to continue? (y/n)" yn
        case $yn in
            [Yy]* ) break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes (y) or no (n).";;
        esac
    done

    rm -r db
    rm -r db_tsv
    rm -r db_temp1
fi

mkdir db

scoref="q"
if [ ! -z "$1" ]
then
    scoref="$1"     #available: q, gg, gg_nhush, q_cc, q_combined
fi

maxconsec=24
if [ ! -z "$2" ]
then
    maxconsec="$2"     # longest consecutive matched stretch
fi

maxid=6
if [ ! -z "$3" ]
then
    maxid="$3"     # max number of identical consecutive nucleotides
fi

targetTemp=72
if [ ! -z "$4" ]
then
    targetTemp="$4"     # minimize melting temperature range around this value
fi

rename 's/.fa$/.hush.out/' HUSH_candidates/*
echo "Building databases."
for f in HUSH_candidates/*; do f2=$(basename $f | sed 's/.hush.out$//'); f3=$(echo "$f2"|sed 's/^sequences_//'); ifpd2 db make -O HUSH_candidates/"$f2".hush.out -T melt/"$f2".tsv -S secs/"$f2".fa.ct db/db."$f3"; done

mkdir db_temp1
echo "Converting to TSV."
for d in db/*; do ifpd2 db dump $d > db_temp1/$(basename $d).tsv; done

mkdir db_temp2
mkdir db_tsv
echo "Attributing oligo score."
echo "Using the following score function: $scoref"
for d in db_temp1/*; do echo $d; sed -r 's/'$'\t''111([0-9]+)987([0-9]+)'$'\t''/'$'\t''\1'$'\t''\2'$'\t''/' $d | sed -r $'s/off_target_no\t/off_target_no\toff_target_sum\t/' > db_temp2/$(basename $d); done

#for d in db_temp1/*; do cat $d | ../escafish_score.py "$scoref" "$2" "$3" > db_tsv/$(basename $d); done
for d in db_temp2/*; do cat $d | ../escafish_score.py "$scoref" "$maxconsec" "$maxid" "$targetTemp" > db_tsv/$(basename $d); done
rm -r db_temp1
rm -r db_temp2

cd ..
