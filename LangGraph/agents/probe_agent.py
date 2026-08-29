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