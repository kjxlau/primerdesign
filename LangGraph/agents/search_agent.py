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
        
        # 2. Group organism terms
        org_query = f'("{organism}"[Organism] OR "{organism}"[All Fields])'
        
        # 3. Exclude WGS (Whole Genome Shotgun) records 
        # (WGS sequences are often unannotated contigs that disrupt alignments)
        exclusions = 'NOT wgs[Property] NOT "whole genome shotgun"[Title]'
        
        # Complete query (Organism AND Gene AND Length NOT WGS)
        # Note: We purposely leave out NOT "partial" to allow partial CDS of target genes
        query = f'{org_query} AND {gene_query} AND 100:15000[SLEN] {exclusions}'
        
        try:
            # Safely fetch IDs using Context Managers
            with Entrez.esearch(db="nucleotide", term=query, retmax=count) as handle:
                search_results = Entrez.read(handle)
                ids = search_results.get("IdList", [])
            
            if not ids: 
                print(f"[SearchAgent] No sequences found for query: {query}")
                return None
            
            print(f"[SearchAgent] Found {len(ids)} sequences. Fetching FASTA data...")
            
            # Fetch the actual FASTA sequence data
            with Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text") as fetch_handle:
                fasta_data = fetch_handle.read()
                return fasta_data
                
        except Exception as e:
            raise Exception(f"NCBI Search Error for query '{query}': {e}")
