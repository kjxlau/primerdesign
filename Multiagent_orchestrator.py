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
        gene_lower = gene.lower()
        if "16s" in gene_lower:
            gene_query = '("16S ribosomal RNA"[Title] OR "16S rRNA"[Title] OR "16S"[Gene])'
        elif "its" in gene_lower or "internal transcribed spacer" in gene_lower:
            gene_query = '("internal transcribed spacer"[Title] OR "ITS"[Title] OR "ITS1"[Title] OR "ITS2"[Title])'
        else:
            gene_query = f'("{gene}"[Gene] OR "{gene}"[Title] OR "{gene}"[All Fields])'
        
        query = f'"{organism}"[Organism] OR "{organism}"[All Fields] OR {gene_query} AND 100:15000[SLEN] NOT "partial"[Title]'
        
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=count)
            ids = Entrez.read(handle)["IdList"]
            if not ids: return None
            fetch_handle = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")
            return fetch_handle.read()
        except Exception as e:
            raise Exception(f"NCBI Search Error: {e}")

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
    def select_probes(self, df_candidates, df_pairs, min_probe_tm, max_results=10, min_spacing=50):
        print(f"[ProbeAgent]: Screening internal probes for {len(df_pairs)} primer pairs...")
        probe_candidates = df_candidates[df_candidates['Tm'] >= min_probe_tm]
        final_sets = []
        df_pairs['Cons_Avg_Pairs'] = (df_pairs['Cons_F'] + df_pairs['Cons_R']) / 2
        sorted_pairs = df_pairs.sort_values(by='Cons_Avg_Pairs', ascending=False)

        for p_idx, pair in sorted_pairs.iterrows():
            too_close = any(abs(pair['Start'] - selected['Start']) < min_spacing for selected in final_sets)
            if too_close: continue

            probe_found = False
            for c_idx, cand in probe_candidates.iterrows():
                if cand['Location'] > (pair['Start'] + pair['AlignWidthF']) and (cand['Location'] + cand['AlignWidth']) < (pair['End'] - pair['AlignWidthR']):
                    final_sets.append({
                        'Forward': pair['Forward'], 'Reverse': pair['Reverse'], 'Probe': cand['Sequence'],
                        'Start': pair['Start'], 'End': pair['End'], 'PLoc': cand['Location'],
                        'TmF': pair['TmF'], 'TmR': pair['TmR'], 'PTm': cand['Tm'], 'Ampsize': pair['Ampsize'], 
                        'Cons_Avg': round((pair['Cons_F'] + pair['Cons_R'] + cand['Conservation'])/3, 1)
                    })
                    probe_found = True
                    break 
            
            if probe_found and len(final_sets) >= max_results:
                break
        return pd.DataFrame(final_sets)

# --- AGENT 6: THE LAMP AGENT ---
class LampAgent:
    """Specializes in isothermal 6-region amplification architecture."""
    def design_lamp(self, df_candidates, min_tm, max_results=10):
        print(f"[LampAgent]: Screening for LAMP primer sets (6 distinct regions)...")
        
        # Sort spatially (Critical for progressive looping!)
        df_filtered = df_candidates[df_candidates['Tm'] >= min_tm].sort_values("Location").reset_index(drop=True)
        records = df_filtered.to_dict('records')
        n = len(records)
        lamp_sets = []

        # Standard Spatial Constraints (LAMP GAP constraints)
        d_F3_F2_min, d_F3_F2_max = 0, 60
        d_F2_F1_min, d_F2_F1_max = 40, 60
        d_F1_B1_min, d_F1_B1_max = 10, 60  # The loop / dumbbell formation region
        d_B1_B2_min, d_B1_B2_max = 40, 60
        d_B2_B3_min, d_B2_B3_max = 0, 60

        for i in range(n):
            r1 = records[i] # F3
            for j in range(i+1, n):
                r2 = records[j] # F2
                dist1 = r2['Location'] - (r1['Location'] + r1['AlignWidth'])
                if dist1 > d_F3_F2_max: break
                if dist1 < d_F3_F2_min: continue
                
                for k in range(j+1, n):
                    r3 = records[k] # F1
                    dist2 = r3['Location'] - (r2['Location'] + r2['AlignWidth'])
                    if dist2 > d_F2_F1_max: break
                    if dist2 < d_F2_F1_min: continue
                    
                    for l in range(k+1, n):
                        r4 = records[l] # B1c
                        dist3 = r4['Location'] - (r3['Location'] + r3['AlignWidth'])
                        if dist3 > d_F1_B1_max: break
                        if dist3 < d_F1_B1_min: continue
                        
                        for m in range(l+1, n):
                            r5 = records[m] # B2c
                            dist4 = r5['Location'] - (r4['Location'] + r4['AlignWidth'])
                            if dist4 > d_B1_B2_max: break
                            if dist4 < d_B1_B2_min: continue
                            
                            for p in range(m+1, n):
                                r6 = records[p] # B3c
                                dist5 = r6['Location'] - (r5['Location'] + r5['AlignWidth'])
                                if dist5 > d_B2_B3_max: break
                                if dist5 < d_B2_B3_min: continue
                                
                                # -- Assemble LAMP Primers --
                                F3_seq = r1['Sequence']
                                F1c = str(Seq(r3['Sequence']).reverse_complement())
                                FIP_seq = F1c + r2['Sequence']
                                B2 = str(Seq(r5['Sequence']).reverse_complement())
                                BIP_seq = r4['Sequence'] + B2
                                B3_seq = str(Seq(r6['Sequence']).reverse_complement())
                                
                                ampsize = (r6['Location'] + r6['AlignWidth']) - r1['Location']
                                avg_cons = round(sum([r['Conservation'] for r in [r1,r2,r3,r4,r5,r6]])/6, 1)
                                
                                lamp_sets.append({
                                    'F3': F3_seq, 'FIP': FIP_seq, 'BIP': BIP_seq, 'B3': B3_seq,
                                    'Start': r1['Location'], 'End': r6['Location'] + r6['AlignWidth'],
                                    'Ampsize': ampsize, 'F3_Tm': r1['Tm'], 'B3_Tm': r6['Tm'],
                                    'Avg_Cons': avg_cons
                                })
                                
                                # Early exit if we have enough sets
                                if len(lamp_sets) >= max_results:
                                    return pd.DataFrame(lamp_sets)
                                    
        return pd.DataFrame(lamp_sets)

# --- THE MASTER ORCHESTRATOR ---
class MasterOrchestrator:
    def __init__(self):
        self.searcher = SearchAgent(os.getenv("NCBI_EMAIL"), os.getenv("NCBI_API_KEY"))
        self.aligner = AlignmentAgent(os.getenv("NCBI_EMAIL"))
        self.analyst = AnalystAgent()
        self.screener = ScreeningAgent()
        self.probe_agent = ProbeAgent()
        self.lamp_agent = LampAgent() # <--- Agent 6 Integrated

    def execute(self, org, gene, count, kw, pcr_params):
        if not isinstance(count, int) or count < 2:
            print(f"[Master]: 'count' of {count} is too low. Defaulting to 10.")
            count = 10

        raw = self.searcher.fetch_sequences(org, gene, count)
        if not raw or raw.count('>') < 2: 
            print(f"[Master]: Error - Not enough sequences found/fetched."); return None

        aln = self.aligner.align(raw)
        if not aln: return None
            
        aln_fn = f"{org}_{gene}_alignment.fasta".replace(" ", "_").lower()
        with open(aln_fn, "w") as f: f.write(aln)
        print(f"[Master]: SUCCESS. Alignment saved to {aln_fn}")
                
        candidates = self.analyst.calculate_stats(aln, kw, 2, 0.05)
        if candidates.empty:
            print("[Master]: No conserved candidate regions found."); return None

        # --- Pipeline A: Standard TaqMan qPCR ---
        pairs = self.screener.screen_pairs(candidates, pcr_params['min_amp'], pcr_params['max_amp'], pcr_params['min_tm'], pcr_params['max_diff'])
        if not pairs.empty:
            final_probe_df = self.probe_agent.select_probes(candidates, pairs, pcr_params['min_probe_tm'])
            if not final_probe_df.empty:
                csv_fn = f"{org}_{gene}_probe_sets.csv".replace(" ","_").lower()
                final_probe_df.to_csv(csv_fn, index=False)
                print(f"[Master]: SUCCESS. TaqMan sets saved to {csv_fn}")
            else: print("[Master]: No valid TaqMan probes found.")
        else: print("[Master]: No valid standard PCR primer pairs found.")

        # --- Pipeline B: LAMP ---
        lamp_df = self.lamp_agent.design_lamp(candidates, pcr_params['min_tm'])
        if not lamp_df.empty:
            lamp_fn = f"{org}_{gene}_lamp_sets.csv".replace(" ","_").lower()
            lamp_df.to_csv(lamp_fn, index=False)
            print(f"[Master]: SUCCESS. LAMP sets saved to {lamp_fn}")
        else:
            print("[Master]: No valid LAMP primer sets found.")

if __name__ == "__main__":
    pcr = {
        'min_amp': int(input("Min Amplicon Size (qPCR): ")), 
        'max_amp': int(input("Max Amplicon Size (qPCR): ")),
        'min_tm': int(input("Min Primer Tm (qPCR/LAMP): ")), 
        'min_probe_tm': int(input("Min Probe Tm (qPCR): ")),
        'max_diff': int(input("Max Length Diff (qPCR): "))
    }
    MasterOrchestrator().execute(
        input("Organism: "), input("Gene: "), int(input("NCBI Records: ")), 
        input("species filter: "), pcr
    )
