#!/bin/bash

cd data
mkdir db

rename 's/.fa$/.hush.out/' HUSH_candidates/*
echo "Building databases"
for f in HUSH_candidates/*; do f2=$(basename $f | sed 's/.hush.out$//'); f3=$(echo "$f2"|sed 's/^sequences_//'); ifpd2 db make -O HUSH_candidates/"$f2".hush.out -T melt/"$f2".tsv -S secs/"$f2".fa.ct db/db."$f3"; done

mkdir db_temp1
echo "Converting to TSV"
for d in db/*; do ifpd2 db dump $d > db_temp1/$(basename $d).tsv; done

mkdir db_temp2
mkdir db_tsv
echo "Attributing oligo score"
for d in db_temp1/*; do echo $d; sed -r 's/'$'\t''111([0-9]+)987([0-9]+)'$'\t''/'$'\t''\1'$'\t''\2'$'\t''/' $d | sed -r $'s/off_target_no\t/off_target_no\toff_target_sum\t/' > db_temp2/$(basename $d); done

for d in db_temp2/*; do cat $d | ../escafish_score.py q > db_tsv/$(basename $d); done
rm -r db_temp1
rm -r db_temp2

cd ..