from typing import Optional
from Bio import Entrez

class SearchAgent:
    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Initializes the NCBI Entrez Search Agent.
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

    def fetch_sequences(self, organism: str, gene: str, count: int = 10) -> Optional[str]:
        """
        Fetches up to `count` fasta sequences from the NCBI Nucleotide database.
        """
        gene_lower = gene.lower()
        
        # 1. Determine Gene Query
        if "16s" in gene_lower:
            gene_query = '("16S ribosomal RNA"[Title] OR "16S rRNA"[Title] OR "16S"[Gene])'
        elif "its" in gene_lower or "internal transcribed spacer" in gene_lower:
            gene_query = '("internal transcribed spacer"[Title] OR "ITS"[Title] OR "ITS1"[Title] OR "ITS2"[Title])'
        else:
            gene_query = f'("{gene}"[Gene] OR "{gene}"[Title] OR "{gene}"[All Fields])'
        
        # 2. Fix Boolean Logic: Group organism terms and use AND
        org_query = f'("{organism}"[Organism] OR "{organism}"[All Fields])'
        
        # Complete query (Organism AND Gene AND Length NOT Partial)
        query = f'{org_query} AND {gene_query} AND 100:15000[SLEN] NOT "partial"[Title]'
        
        try:
            # 3. Use Context Managers (with statement) to ensure handles auto-close
            with Entrez.esearch(db="nucleotide", term=query, retmax=count) as handle:
                search_results = Entrez.read(handle)
                # Safely get IdList, defaulting to an empty list if missing
                ids = search_results.get("IdList", [])
            
            if not ids: 
                print(f"No sequences found for query: {query}")
                return None
            
            # Fetch the actual sequence data
            with Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text") as fetch_handle:
                fasta_data = fetch_handle.read()
                return fasta_data
                
        except Exception as e:
            # Better error message reporting exactly what query failed
            raise Exception(f"NCBI Search Error for query '{query}': {e}")