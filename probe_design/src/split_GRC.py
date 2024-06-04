from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os

'''
Split a GRC genome fasta file into chromosomes.
The chromosome id is extracted from the fasta title.
The output files are named as "save_loc/prefix.chromosome.chr_id.fa"
Also, FASTA sequences are transformed to uppercase.
'''
all_chr: list[int|str] = list(range(1,23))+['X','Y','MT'] # list of human chromosomes
all_chr = [str(c) for c in all_chr]

def split_fasta(genome_fasta:str, prefix:str|None=None,save_loc:str='.')->None:
    '''
    Read a genome fasta file and split it into chromosomes.
    The chromosome id is extracted from the fasta title.
    The output files are named as "save_loc.prefix.chromosome.chr_id.fa"
    '''

    if prefix is None:
        position = genome_fasta.rfind('.primary_assembly')
        prefix = genome_fasta[:position]
  
    with open(genome_fasta, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            title = record.description
            sequence = record.seq.upper() # Convert to uppercase
            words = title.split(' ')

            chr_id = words[0] # first word is the chromosome id
            if chr_id not in all_chr:
                continue

            fname = os.path.join(save_loc, f'{prefix}.chromosome.{chr_id}.fa') # Output file name
            with open(fname, 'w') as out:
                # Create a new SeqRecord for the sequence
                print(f'\033[92mWriting chromosome {chr_id} to {fname}\033[0m',end=' ')
                new_record = SeqRecord(sequence, id=record.id, description=title)
                # Write the SeqRecord to the file
                SeqIO.write(new_record, out, 'fasta')   # 60 characters per line
                print('\u2705')
        
    return

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Split a GRC genome fasta file into chromosomes.')
    argparser.add_argument('-g','--genome', help='Input genome fasta file.',default='data/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa')
    argparser.add_argument('-p','--prefix', help='Prefix for the output files.', default=None)
    argparser.add_argument('-s','--save_loc', help='Location to save the output chromosome files.', default='.')
    args = argparser.parse_args()
    split_fasta(args.genome, args.prefix, args.save_loc)
