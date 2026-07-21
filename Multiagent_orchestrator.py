import os
import time
import requests
import pandas as pd
from io import StringIO
from dotenv import load_dotenv
from Bio import Entrez, AlignIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Align import AlignInfo, MultipleSeqAlignment

load_dotenv()

# --- AGENT 1: THE SEARCH AGENT ---
class SearchAgent:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key

    def fetch_sequences(self, organism, gene, count):
        # 1. Expand standard abbreviations if it's 16S
        if "16s rrna" in gene.lower():
            gene_query = '("16S ribosomal RNA"[Title] OR "16S rRNA"[Title] OR "16S"[Gene])'
        else:
            # For standard genes (e.g. gyrB, recA), check both Gene and Title annotations
            gene_query = f'("{gene}"[Gene] OR "{gene}"[Title])'

        # 2. Build the robust query
        # [SLEN] strictly limits results to 100bp - 15,000bp, completely eliminating whole genomes (millions of bp)
        # without needing fragile "NOT genome" text matching.
        query = f'"{organism}"[Organism] AND {gene_query} AND 100:15000[SLEN] NOT "partial"[Title]'
        
        print(f"\n[SearchAgent]: Querying NCBI for: {query}")
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=count)
            ids = Entrez.read(handle)["IdList"]
            if not ids: 
                return None
            fetch_handle = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")
            return fetch_handle.read()
        except Exception as e:
            print(f"[SearchAgent]: Error - {e}")
            return None
            
# --- AGENT 2: THE ALIGNMENT AGENT ---
class AlignmentAgent:
    def __init__(self, email):
        self.email = email
        self.url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"

    def align(self, fasta_content):
        print("[AlignmentAgent]: Submitting to EBI Clustal Omega...")
        try:
            job_res = requests.post(f"{self.url}/run", data={'email': self.email, 'sequence': fasta_content, 'stype': 'dna', 'outfmt': 'fa'})
            if not job_res.ok:
                print(f"[AlignmentAgent]: API Error - {job_res.text}")
                return None
            job_id = job_res.text
            while True:
                status = requests.get(f"{self.url}/status/{job_id}").text
                if status == "FINISHED": break
                elif status in ["FAILURE", "NOT_FOUND"]: return None
                print(f"[AlignmentAgent]: Status: {status}...", end="\r")
                time.sleep(5)
            return requests.get(f"{self.url}/result/{job_id}/aln-fasta").text
        except Exception as e:
            print(f"[AlignmentAgent]: Error - {e}"); return None

# --- AGENT 3: THE ANALYST AGENT ---
class AnalystAgent:
    def calculate_stats(self, alignment_text, target_kw, degen_limit, cutoff):
        print(f"[AnalystAgent]: Finding conserved candidates for '{target_kw}'...")
        handle = StringIO(alignment_text)
        alignment = AlignIO.read(handle, "fasta")
        target_list = [r for r in alignment if target_kw.lower() in r.description.lower()]
        target_msa = MultipleSeqAlignment(target_list) if target_list else alignment
        summary = AlignInfo.SummaryInfo(target_msa)
        consensus = str(summary.dumb_consensus(threshold=0.7, ambiguous='N'))
        results = []
        for width in range(20, 31):
            for i in range(len(consensus) - width + 1):
                window = consensus[i:i+width].replace('-', '')
                if len(window) < 18 or window.count('N') > degen_limit: continue
                matches = sum(1 for r in target_msa if str(r.seq[i:i+width]).replace('-', '') == window)
                ratio = matches / len(target_msa)
                if ratio >= (1 - cutoff):
                    results.append({
                        'Sequence': window, 'Location': i, 'Length': len(window), 'AlignWidth': width, 
                        'GC': round((window.count('G') + window.count('C')) / len(window) * 100, 1),
                        'Tm': round(mt.Tm_NN(Seq(window)), 1), 'Conservation': round(ratio * 100, 1)
                    })
        return pd.DataFrame(results)

# --- AGENT 4: THE SCREENING AGENT (Primer Pairs) ---
class ScreeningAgent:
    def screen_pairs(self, df_candidates, min_amp, max_amp, min_tm, max_len_diff):
        print(f"[ScreeningAgent]: Pairing Forward and Reverse primers...")
        valid_indices = df_candidates[df_candidates['Tm'] >= min_tm].index.tolist()
        pairs = []
        for i in valid_indices:
            for j in valid_indices:
                fwd, rev = df_candidates.loc[i], df_candidates.loc[j]
                amplicon_size = (rev['Location'] + rev['Length']) - fwd['Location']
                if min_amp <= amplicon_size <= max_amp and abs(fwd['Length'] - rev['Length']) <= max_len_diff:
                    if rev['Location'] > (fwd['Location'] + fwd['AlignWidth']):
                        pairs.append({
                            'Forward': fwd['Sequence'], 'Reverse': str(Seq(rev['Sequence']).reverse_complement()),
                            'AlignWidthF': fwd['AlignWidth'], 'AlignWidthR': rev['AlignWidth'],
                            'Start': fwd['Location'], 'End': rev['Location'] + rev['Length'],
                            'LenF': fwd['Length'], 'LenR': rev['Length'], 'Ampsize': amplicon_size,
                            'TmF': fwd['Tm'], 'TmR': rev['Tm'], 'GCF': fwd['GC'], 'GCR': rev['GC'],
                            'Cons_F': fwd['Conservation'], 'Cons_R': rev['Conservation']
                        })
        return pd.DataFrame(pairs)

# --- AGENT 5: THE PROBE AGENT (qPCR Specialist) ---
class ProbeAgent:
    """Specializes in finding internal probes between primer pairs."""
    def select_probes(self, df_candidates, df_pairs, min_probe_tm, max_results=10):
        print(f"[ProbeAgent]: Screening internal probes for {len(df_pairs)} primer pairs...")
        # 1. Get candidates that meet probe Tm requirements
        probe_candidates = df_candidates[df_candidates['Tm'] >= min_probe_tm]
        final_sets = []
        count = 0

        for p_idx, pair in df_pairs.iterrows():
            for c_idx, cand in probe_candidates.iterrows():
                # Logic: Probe must start after Fwd ends, and end before Rev starts
                probe_start = cand['Location']
                probe_end = cand['Location'] + cand['AlignWidth']
                fwd_end = pair['Start'] + pair['AlignWidthF']
                rev_start = pair['End'] - pair['AlignWidthR']

                if probe_start > fwd_end and probe_end < rev_start:
                    count += 1
                    final_sets.append({
                        'Forward': pair['Forward'], 'Reverse': pair['Reverse'], 'Probe': cand['Sequence'],
                        'Start': pair['Start'], 'End': pair['End'], 'PLoc': cand['Location'],
                        'TmF': pair['TmF'], 'TmR': pair['TmR'], 'PTm': cand['Tm'],
                        'Ampsize': pair['Ampsize'], 'Conservation_Avg': (pair['Cons_F'] + pair['Cons_R'] + cand['Conservation'])/3
                    })
                    if count >= max_results: break
            if count >= max_results: break
            
        return pd.DataFrame(final_sets)

# --- THE MASTER ORCHESTRATOR ---
class MasterOrchestrator:
    def __init__(self):
        self.searcher = SearchAgent(os.getenv("NCBI_EMAIL"), os.getenv("NCBI_API_KEY"))
        self.aligner = AlignmentAgent(os.getenv("NCBI_EMAIL"))
        self.analyst = AnalystAgent()
        self.screener = ScreeningAgent()
        self.probe_agent = ProbeAgent()

    def execute(self, org, gene, count, kw, pcr_params):
        # 1. Enforce a minimum count for Multiple Sequence Alignment
        if not isinstance(count, int) or count < 2:
            print(f"[Master]: 'count' of {count} is too low for alignment. Defaulting to 10.")
            count = 10

        raw = self.searcher.fetch_sequences(org, gene, count)
        if not raw: 
            return None
            
        # 2. Validate that NCBI actually returned multiple FASTA sequences
        if raw.count('>') < 2:
            print(f"[Master]: Error - NCBI returned less than 2 sequences. Cannot perform alignment.")
            return None

        # 3. Perform Alignment
        aln = self.aligner.align(raw)
        if not aln: 
            return None
            
        # Save the alignment to a FASTA file 
        aln_fn = f"{org}_{gene}_alignment.fasta".replace(" ", "_").lower()
        with open(aln_fn, "w") as f:
            f.write(aln)
        print(f"[Master]: SUCCESS. Alignment saved to {aln_fn}")
                
        # 4. Analyze candidates
        candidates = self.analyst.calculate_stats(aln, kw, 2, 0.05)
        
        if candidates.empty:
            print("[Master]: No conserved candidate regions found. Try relaxing the conservation cutoff.")
            return None

        # 5. Screen pairs
        pairs = self.screener.screen_pairs(candidates, pcr_params['min_amp'], pcr_params['max_amp'], pcr_params['min_tm'], pcr_params['max_diff'])
        
        if pairs.empty:
            print("[Master]: No primer pairs found."); return None

        # 6. Select Probes
        final_df = self.probe_agent.select_probes(candidates, pairs, pcr_params['min_probe_tm'])

        # 7. Output Final CSV
        if not final_df.empty:
            csv_fn = f"{org}_{gene}_probe_sets.csv".replace(" ","_").lower()
            final_df.to_csv(csv_fn, index=False)
            print(f"[Master]: SUCCESS. TaqMan sets saved to {csv_fn}")
            print(final_df.head(5))
            return csv_fn
        else:
            print("[Master]: No valid probes found between the primer pairs.")
            return None

if __name__ == "__main__":
    pcr = {
        'min_amp': int(input("Min Amplicon Size: ")), 'max_amp': int(input("Max Amplicon Size: ")),
        'min_tm': int(input("Min Primer Tm: ")), 'min_probe_tm': int(input("Min Probe Tm: ")),
        'max_diff': int(input("Max Length Diff: "))
    }
    MasterOrchestrator().execute(
        input("Organism: "), input("Gene: "), int(input("NCBI Records: ")), 
        input("species: "), pcr
    )
