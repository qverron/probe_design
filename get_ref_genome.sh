 # Download ref genomes
prefix="Homo_sapiens.GRCh37.dna"
address="http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/"
total=22
# ex: get-ref-genome.sh 

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

