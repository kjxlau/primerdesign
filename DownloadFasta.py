from Bio import Entrez

Entrez.email = "kennyjxlau@gmail.com"
Entrez.api_key = "6b675e217b0d870e6ca896cd257b4c073708"

no_of_records=input("Enter No of Records to download: ")
organism=input("Enter organism of interest: ")
#DENV1: txid11053[Organism:exp]
#DENV2: txid11060[Organism:exp]
#DENV3: txid11069[Organism:exp]
#DENV4: txid11070[Organism:exp]
gene=input("Enter target gene: ")
country=input("Country of interest: ")
search_term=organism+str(" [organism] AND ")+gene+str(" [gene] AND ")+country+str(" [All Fields]")

IDs = Entrez.read(Entrez.esearch(db="nucleotide", retmax=no_of_records, term=search_term)) ["IdList"]
filename=input("Filename to save as fasta: ")+str(".fasta")
with open(filename, 'w') as f:
    for ID in IDs:
        seq=Entrez.efetch(db="nucleotide", id=ID, rettype="fasta", retmode="text").read()
        print(seq)
        f.write(seq)
