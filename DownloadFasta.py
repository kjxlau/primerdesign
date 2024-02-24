from Bio import Entrez

Entrez.email = "kennyjxlau@gmail.com"
Entrez.api_key = "6b675e217b0d870e6ca896cd257b4c073708"

IDs = Entrez.read(Entrez.esearch(db="nucleotide", retmax=10, term="Dengue AND Singapore AND China AND India AND Indonesia")) ["IdList"]
with open('DENV_Asia.fasta', 'w') as f:
    for ID in IDs:
        seq=Entrez.efetch(db="nucleotide", id=ID, rettype="fasta", retmode="text").read()
        print(seq)
        f.write(seq)

