from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from tqdm import tqdm


argparser = argparse.ArgumentParser(description="Split T2T reference into individual chromosomes")
argparser.add_argument("-g",'--genome', help="Path to the reference fasta file", type=str, default="~/ref/T2T.fa")
argparser.add_argument('-n','--name',help="Name of the reference", type = str, default="CHM13.T2T")
argparser.add_argument('-d','--directory',help="Path to the output directory", type = str, default="./data/ref")
args = argparser.parse_args()

def extract_headers_and_sequences(fasta_file:str,ref_name:str = "CHM13.T2T", dir=None) -> tuple[list[str], list[str], list[str]]:
    headers = []
    sequences = []

    if dir:
        ref_name = f"{dir}/{ref_name}"

    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(record.description)
        sequences.append(record.seq)

    chr_ids = [x.split(" ")[6].strip(',') for x in headers]
    file_names = [f"{ref_name}.chromosome.{x}.fa" for x in chr_ids]

    return file_names, headers, sequences


def save_fasta_files(file_names:list[str], headers:list[str], sequences:list[str]) -> None:


    for i in tqdm(range(len(file_names))):
        sequence = Seq(sequences[i])
        record = SeqRecord(sequence, id=headers[i])
        with open(file_names[i], "w") as f:
            SeqIO.write(record, f, "fasta")
    return


if __name__ == "__main__":
    ref = 'GCF_009914755.1_T2T-CHM13v2.0_genomic.fna'
    save_fasta_files(
        *extract_headers_and_sequences(
            fasta_file=args.genome,
            dir=args.directory,
            ref_name=args.name,
            )
        )



