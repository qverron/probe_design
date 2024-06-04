import requests
from bs4 import BeautifulSoup



# get the website
ws = requests.get("https://www.ensembl.org/{Homo_sapiens}/Info/Index".format())

soup = BeautifulSoup(ws.content, 'html.parser')
title = soup.find('title').string

# get the species and build number
species = title.split()[0]
build = title.split()[-1]
