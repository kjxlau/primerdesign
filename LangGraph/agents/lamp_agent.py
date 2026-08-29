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
