
# Python code to download the chromosome 17 of the human genome from https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz

import requests
import os
import gzip

# Human (homo_sapiens) and Mouse (mus_musculus) genomes respectively
GRCh = "https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh{build}.dna.chromosome.{chr}.fa.gz"
GRCm = "https://ftp.ensembl.org/pub/release-{release}/fasta/mus_musculus/dna/Mus_musculus.GRCm{build}.dna.chromosome.{chr}.fa.gz"

# available only mouse and human
available = ['human','homo_sapiens','mouse','mus_musculus']
human = ['human','homo_sapiens']
mouse = ['mouse','mus_musculus']

human_all_chr: list[int|str] = list(range(1,23))+['X','Y'] # list of human chromosomes
mouse_all_chr: list[int|str] = list(range(1,19))+['X','Y'] # list of mouse chromosomes

def download_chr(chr_folder:os.PathLike="data/ref",
                 chr:int|str=17,
                 release:int|str=109,
                 build:int=38,
                 species:str="human",
                 unzip:bool=True,
                 remove_gz:bool=False,
                 file_name:str|None=None,
                 url_chr:str|None=None)->None:
    
    """
    Downloads the specified chromosome of a genome from Ensembl FTP server. Requires internet connection...

    Args:
        chr_folder (os.PathLike, optional): 
            Path to the folder where the chromosome file will be saved. Defaults to "data/ref".
        chr (int|str, optional): 
            Chromosome number or identifier. Defaults to 17.
        release (int|str, optional): 
            Release version of the genome. Defaults to 109.
        build (int, optional):
            Build version of the genome. Defaults to 38.
        species (str, optional):
            Species of the genome (human or mouse). Defaults to "human".
        unzip (bool, optional): 
            Flag to indicate whether to unzip the downloaded file. Defaults to True.
        remove_gz (bool, optional): 
            Flag to indicate whether to remove the downloaded gzipped file. Defaults to False.

    Returns:
        None
    """

    # matching species with the link
    species = species.lower()
    if url_chr is None:
        if species in human:
            url_chr = GRCh
        elif species in mouse:
            url_chr = GRCm
        else:
            print(f"species should be one of {available}")
            print(f"You entered {species}")
            return
        print(f"downloading {species} chromosome {chr}...")
        url_chr = url_chr.replace('{release}',f'{release}'
                                    ).replace('{chr}',f'{chr}'
                                        ).replace('{build}',f'{build}') # placing release and chr numbers
    else:
        print(f"downloading from {url_chr}...")
        print(f"if reference genome, this may take some time: 5-15 mins depending on your internet speed")
        pass # url_chr = url_chr

    if file_name is None:
        file_name = url_chr.split('/')[-1] # last part of the link will be the file name
    else:
        pass # file_name = file_name
    
    # creates chr_folder if not exists
    if not os.path.exists(chr_folder):
        try: # just in case the parent directories do not exist
            os.mkdir(chr_folder)
        except:
            print(f'parent directory {os.path.dirname(chr_folder)} does not exist')
            print('exiting...')
            return
    
    chr_path = os.path.join(chr_folder,file_name)
    
    # read the link to a variable chr_fasta_gz
    if not os.path.exists(chr_path): # do not overwrite if it exists
        
        try:
            chr_fasta_gz = requests.get(url_chr)
            open(chr_path,"wb").write(chr_fasta_gz.content)# write chr_fasta_gz into the chromosome path
        except:
            print("Link does not work:",url_chr)
            print("Exiting...")
            return
    else:
        print(f'file {chr_path} already exists. Skipping...')


    if unzip: # unzip the .gz file if requested
        chr_fasta = chr_path[:-3]  # Remove the '.gz' extension
        if not os.path.exists(chr_fasta): # do not overwrite chromosome fasta file if exists
            print(f'unzipping...')
            # Unzip .gz file and save the unzipped content to a new file
            open(chr_fasta, 'wb').write(gzip.open(chr_path, 'rb').read())
            print(f'unzipped file saved to {chr_fasta}')
        else:
            print(f'file {chr_fasta} already exists. Skipping...')

    if remove_gz: # remove the .gz file if requested
        if os.path.exists(chr_path):
            print(f'Removing gzipped file {chr_path} ...')
            os.remove(chr_path)
            print('Done')
        else: # could already be removed
            print(f'file {chr_path} does not exist. Skipping...')


    return

def download_chr_list(list_chr:str|list = 'all',   
                      chr_folder:os.PathLike="data/ref",
                      release:int|str=109,
                      build:int=38,
                      species:str="human",
                      unzip:bool=True,
                      remove_gz:bool=False)->None:
    
    """
    Downloads the specified chromosome lists of a genome from Ensembl FTP server. 
        Requires internet connection...

    Args:
        list_chr (list[int|str], optional): 
            List of chromosome numbers or identifiers. Defaults to list(range(1,23))+['X','Y'].
        chr_folder (os.PathLike, optional): 
            Path to the folder where the chromosome file will be saved. Defaults to "data/ref".
        release (int|str, optional): 
            Release version of the genome. Defaults to 109.
        build (int, optional):
            Build version of the genome. Defaults to 38.
        species (str, optional):
            Species of the genome (human or mouse). Defaults to "human".
        unzip (bool, optional): 
            Flag to indicate whether to unzip the downloaded file. Defaults to True.
        remove_gz (bool, optional): 
            Flag to indicate whether to remove the downloaded gzipped file. Defaults to False.

    Returns:
        None
    """
    if list_chr=="all":
        if species in human :
            list_chr = human_all_chr
        elif species in mouse:
            list_chr = mouse_all_chr
        else:
            print(f"species should be one of {available}")
            print(f"You entered {species}")
            return
    elif type(list_chr) == list or type(list_chr) == tuple:
        pass # 
    else:
        print('list_chr should either be "all", a list or tuple')
        print(f'You entered {type(list_chr)} type which is not comptatible..')
        return

    for chr in list_chr:
        print(f'downloanding... {species=} {chr=} {release=} {build=} ')
        download_chr(chr_folder=chr_folder,
                     chr=chr,
                     release=release,
                     species=species,
                     unzip=unzip,
                     remove_gz=remove_gz,
                     )
    return

def download_ref_genome(build:int=38,
                        release:int=109,
                        species:str="human",
                        chr_folder:os.PathLike="data/ref",
                        file_name:str='genome.fa.gz',
                        unzip:bool=True,
                        remove_gz:bool=False)->None:
    
    if species in human:
        url_chr = GRCh.replace('{release}',f'{release}').replace('{build}',f'{build}')
        
    elif species in mouse:
        url_chr = GRCm.replace('{release}',f'{release}').replace('{build}',f'{build}')


    url_chr = url_chr.replace('chromosome.{chr}','primary_assembly')
    download_chr(chr_folder=chr_folder,
                 url_chr=url_chr,
                 unzip=unzip,
                 file_name=file_name,
                 remove_gz=remove_gz)
    return



if __name__ == "__main__":
    download_ref_genome() # works for the default values