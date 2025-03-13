from Bio import Entrez

Entrez.email = "kennyjxlau@gmail.com"
Entrez.api_key = "6b675e217b0d870e6ca896cd257b4c073708"

no_of_records=input("Enter No of Records to download: ")
organism=input("Enter organism of interest (search strings): ")

search_term=organism
print(search_term)

IDs = Entrez.read(Entrez.esearch(db="nucleotide", retmax=no_of_records, term=search_term)) ["IdList"]
filename=input("Filename to save as fasta: ")+str(".fasta")
with open(filename, 'w') as f:
    for ID in IDs:
        seq=Entrez.efetch(db="nucleotide", id=ID, rettype="fasta", retmode="text").read()
        print(seq)
        f.write(seq)
