
# Python code to download the chromosomes and genome of a species from Ensembl FTP server
import requests
import os
import gzip
import argparse
from bs4 import BeautifulSoup
from .split_GRC import split_fasta
from tqdm import tqdm

# Human (homo_sapiens) and Mouse (mus_musculus) genomes respectively
GRCh = "https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh{build}.dna.chromosome.{chr}.fa.gz"
GRCm = "https://ftp.ensembl.org/pub/release-{release}/fasta/mus_musculus/dna/Mus_musculus.GRCm{build}.dna.chromosome.{chr}.fa.gz"

# available only mouse and human
available = ['homo_sapiens','mus_musculus']
human = ['homo_sapiens']
mouse = ['mus_musculus']

human_all_chr: list[int|str] = list(range(1,23))+['X','Y','MT'] # list of human chromosomes
mouse_all_chr: list[int|str] = list(range(1,19))+['X','Y','MT'] # list of mouse chromosomes

#args
argparser = argparse.ArgumentParser(prog='download ensemble genome',description="download ensemble genome")
argparser.add_argument('-s','--species',type=str,default='homo_sapiens',choices=["homo_sapiens","mus_musculus"])
argparser.add_argument('-b','--build',type=int,default=38,help="the build number of the genome")
argparser.add_argument('-r','--release',default=109,help="release number of the build")
argparser.add_argument('-f','--filename',default=None,help="give a specific name to the downloaded file ")
argparser.add_argument('-k','-keep',action='store_true',help="whether to keep gzip files")
argparser.add_argument('-split',action='store_true',help='whether to split into chromosomes')
args = argparser.parse_args()


def get_ensembl_release(species:str='Homo_sapiens')->int:
    '''
    get the latest release number for ensemble
    '''

    # get the website
    ws = requests.get(f"https://www.ensembl.org/{species}/Info/Index")

    # get the tile
    soup = BeautifulSoup(ws.content, 'html.parser')
    title = soup.find('title').string

    # get the species and build number
    species = title.split()[0]
    latest_release = title.split()[-1]

    print(f"\033[1m\033[94m build (latest) {latest_release} for {args.species} \033[0m")

    return latest_release

def download_chr(chr_folder:os.PathLike="data/ref",
                 chr:int|str=17,
                 release:int|str=109,
                 build:int=38,
                 species:str="homo_sapiens",
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
            Species of the genome (homo_sapiens or mus_musculus). Defaults to "homo_sapiens".
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
        print(f"\033[91m downloading from\033[0m\n{url_chr} ...")
        url_chr = url_chr.replace('{release}',f'{release}'
                                    ).replace('{chr}',f'{chr}'
                                        ).replace('{build}',f'{build}') # placing release and chr numbers
    else:
        print(f"downloading from {url_chr} ...")
        if 'assembly' in url_chr:
            print(f"\033[91mREFERENCE genome, this may take some time: 5-15 mins depending on your internet speed\033[0m")
        pass # url_chr = url_chr

    if file_name is None:
        file_name = url_chr.split('/')[-1] # last part of the link will be the file name
    else:
        pass # file_name = file_name
    
    # creates chr_folder if not exists
    if not os.path.exists(chr_folder):
        try: # just in case the directories do not exist
            os.makedirs(chr_folder,exist_ok=True)
        except:
            print(f'parent directory {os.path.dirname(chr_folder)} does not exist')
            print('exiting...')
            return
    
    chr_path = os.path.join(chr_folder,file_name)
    
    # read the link to a variable chr_fasta_gz
    if not os.path.exists(chr_path): # do not overwrite if it exists
        
        #try:
        if True:
            chr_fasta_gz = requests.get(url_chr,stream=True)
            total_size = int(chr_fasta_gz.headers.get('content-length', 0))
            block_size = 1024 #1 Kbyte
            with tqdm(total=total_size, unit='B', unit_scale=True) as progress_bar:
                with open(chr_path,"wb") as f:
                    for data in chr_fasta_gz.iter_content(block_size):
                        progress_bar.update(len(data))
                        f.write(data)
            #open(chr_path,"wb").write(chr_fasta_gz.content)# write chr_fasta_gz into the chromosome path
        #except:
        #    print("Link does not work:",url_chr)
        #    print("Exiting...")
         #   return
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
            print(f'Removing gzipped file {chr_path} ...', end=' ')
            os.remove(chr_path)
            print('\u2705')
        else: # could already be removed
            print(f'file {chr_path} does not exist. Skipping...')


    return

def download_chr_list(list_chr:str|list = 'all',   
                      chr_folder:os.PathLike="data/ref",
                      release:int|str=109,
                      build:int=38,
                      species:str="homo_sapiens",
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
            Species of the genome (homo_sapiens or mus_musculus). Defaults to "homo_sapiens".
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
                        species:str="homo_sapiens",
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

def main():

    if args.release != 'latest':
        download_ref_genome(
            build=args.build,
            release=args.release,
            species=args.species,
            chr_folder="data/ref",
            file_name=args.filename,
            unzip=True,
            remove_gz=not args.k
        )
    else:
        release = get_ensembl_release(species=args.species.capitalize())
        download_ref_genome(
            build=args.build,
            release=release,
            species=args.species,
            chr_folder="data/ref",
            file_name=args.filename,
            unzip=True,
            remove_gz=not args.k
        )

    if args.split:
        import glob
        fasta = glob.glob('data/ref/*.dna.primary_assembly.fa')[0]
        split_fasta(genome_fasta=fasta,prefix=None,save_loc='.') # split the genome into chromosomes
    

if __name__ == "__main__":
    main()

