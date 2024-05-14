from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os

'''
Split a CHM13 T2T (telomere-to-telomere) genome fasta file into chromosomes.
The chromosome id is extracted from the fasta title.
The output files are named as "save_loc/prefix.chromosome.chr_id.fa"
Also, FASTA sequences are transformed to uppercase.
'''

argparser = argparse.ArgumentParser(description='Split a CHM13 T2T (telomere-to-telomere) genome fasta file into chromosomes.')
argparser.add_argument('-g','--genome', help='Input genome fasta file.',default='GCA_009914755.4_T2T-CHM13v2.0_genomic.fna')
argparser.add_argument('-p','--prefix', help='Prefix for the output files.', default='CHM13.T2T')
argparser.add_argument('-s','--save_loc', help='Location to save the output chromosome files.', default='.')
args = argparser.parse_args()

def split_fasta(genome_fasta:str, prefix:str='chm13',save_loc:str='.')->None:
    '''
    Read a genome fasta file and split it into chromosomes.
    The chromosome id is extracted from the fasta title.
    The output files are named as "save_loc.prefix.chromosome.chr_id.fa"
    '''
    with open(genome_fasta, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            title = record.description
            sequence = record.seq.upper() # Convert to uppercase
            words = title.split(' ')
            if not 'mitochondrion' in title: # Not a mitochondrial chromosome
                chr_id = words[-1] # last word is the chromosome id
            else : # Mitochondrial chromosome
                chr_id = 'M' 

            fname = os.path.join(save_loc, f'{prefix}.chromosome.{chr_id}.fa') # Output file name
            with open(fname, 'w') as out:
                # Create a new SeqRecord for the sequence
                print(f'Writing chromosome {chr_id} to {save_loc}{prefix}.chromosome.{chr_id}.fa ...',end='  ')
                new_record = SeqRecord(sequence, id=record.id, description=title)
                # Write the SeqRecord to the file
                SeqIO.write(new_record, out, 'fasta')   # 60 characters per line
                print(f' Done!')


if __name__ == '__main__':
    split_fasta(args.genome, args.prefix, args.save_loc)
    test_script = 'python3 split_it.py -g GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -p chm13 -s ./'
    test_for_hpc = 'python3 split_it.py -g data/ref/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -p CHM13.T2T -s data/ref/'