#!/bin/bash

mkdir -p data/ref
cd data/ref
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
gunzip GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
prb split_T2T -g GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -p CHM13.T2T -s .
mv GCA_009914755.4_T2T-CHM13v2.0_genomic.fna genome.fa