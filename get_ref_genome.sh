#!/bin/bash

# Download ref genomes
prefix="Mus_musculus.GRCm38.dna"
address="http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/"
total=19

# Human GRCh37:
#prefix="Homo_sapiens.GRCh37.dna"
#address="https://ftp.ensembl.org/pub/grch37/release-109/fasta/homo_sapiens/dna/"
#total=22

# Human GRCh38:
#prefix="Homo_sapiens.GRCh38.dna"
#address="https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/“
#total=22

mkdir -p data/ref
cd data
echo $(seq 1 "$total") X Y MT | tr " " "\n" | xargs -I% echo chromosome.% > ref/chrs.list
echo nonchromosomal >> ref/chrs.list

for chr in $(cat ref/chrs.list)
do
     wget -O ref/"$prefix".${chr}.noChr.fa.gz "$address"/"$prefix".${chr}.fa.gz
     zcat ref/"$prefix"."$chr".noChr.fa.gz | sed 's/>/>chr/' | cut -d" " -f 1 > ref/"$prefix".${chr}.fa
     rm ref/"$prefix".${chr}.noChr.fa.gz    
done

