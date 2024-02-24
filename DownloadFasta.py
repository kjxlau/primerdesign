from Bio import Entrez

Entrez.email = "kennyjxlau@gmail.com"
Entrez.api_key = "6b675e217b0d870e6ca896cd257b4c073708"

IDs = Entrez.read(Entrez.esearch(db="nucleotide", retmax=10, term="txid11053[Organism:exp] AND Thailand[All Fields]")) ["IdList"]
#DENV1: txid11053[Organism:exp]
#DENV2: txid11060[Organism:exp]
#DENV3: txid11069[Organism:exp]
#DENV4: txid11070[Organism:exp]
with open('DENV1_Thailand.fasta', 'w') as f:
    for ID in IDs:
        seq=Entrez.efetch(db="nucleotide", id=ID, rettype="fasta", retmode="text").read()
        print(seq)
        f.write(seq)
