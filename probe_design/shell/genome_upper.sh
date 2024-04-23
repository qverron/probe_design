#!/bin/bash


# Convert the reference genome to uppercase for compatibility
cd data/ref
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' genome.fa > genome.fa.2
rm genome.fa
mv genome.fa.2 genome.fa