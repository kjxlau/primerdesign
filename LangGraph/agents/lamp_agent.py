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
    """Specializes in isothermal 8-region amplification architecture (Includes Loop Primers)."""
    
    def design_lamp(self, df_candidates, min_tm, max_results=10):
        print(f"[LampAgent]: Screening for complete LAMP primer sets (6 primers, 8 distinct regions)...")
        
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
                    
                    # --- NEW: Identify LoopF (LF) candidates in the gap between F2 and F1 ---
                    valid_LFs = []
                    for x in range(j+1, k):
                        r_LF = records[x]
                        # LF must fit completely within the gap between F2's 3' end and F1's 5' end
                        if (r_LF['Location'] >= r2['Location'] + r2['AlignWidth']) and \
                           (r_LF['Location'] + r_LF['AlignWidth'] <= r3['Location']):
                            valid_LFs.append(r_LF)
                    if not valid_LFs: 
                        continue # Skip if no LoopF fits (enforce full 6-primer set)
                    
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
                            
                            # --- NEW: Identify LoopB (LB) candidates in the gap between B1c and B2c ---
                            valid_LBs = []
                            for y in range(l+1, m):
                                r_LB = records[y]
                                # LB must fit completely within the gap between B1c's 3' end and B2c's 5' end
                                if (r_LB['Location'] >= r4['Location'] + r4['AlignWidth']) and \
                                   (r_LB['Location'] + r_LB['AlignWidth'] <= r5['Location']):
                                    valid_LBs.append(r_LB)
                            if not valid_LBs:
                                continue # Skip if no LoopB fits
                            
                            for p in range(m+1, n):
                                r6 = records[p] # B3c
                                dist5 = r6['Location'] - (r5['Location'] + r5['AlignWidth'])
                                if dist5 > d_B2_B3_max: break
                                if dist5 < d_B2_B3_min: continue
                                
                                # -- Assemble Complete Accelerated LAMP Primers --
                                best_LF = valid_LFs[0] # Select the first viable LoopF
                                best_LB = valid_LBs[0] # Select the first viable LoopB
                                
                                F3_seq = r1['Sequence']
                                
                                # FIP = F1c + F2
                                F1c = str(Seq(r3['Sequence']).reverse_complement())
                                FIP_seq = F1c + r2['Sequence']
                                
                                # LoopF targets the loop on the sense strand, so it must be Antisense
                                LF_seq = str(Seq(best_LF['Sequence']).reverse_complement())
                                
                                # LoopB targets the loop on the antisense strand, so it must be Sense
                                LB_seq = best_LB['Sequence']
                                
                                # BIP = B1c + B2
                                B2 = str(Seq(r5['Sequence']).reverse_complement())
                                BIP_seq = r4['Sequence'] + B2
                                
                                # B3 = Reverse Complement of B3c
                                B3_seq = str(Seq(r6['Sequence']).reverse_complement())
                                
                                ampsize = (r6['Location'] + r6['AlignWidth']) - r1['Location']
                                
                                # Average conservation score across all 8 regions
                                regions = [r1, r2, best_LF, r3, r4, best_LB, r5, r6]
                                avg_cons = round(sum([r['Conservation'] for r in regions]) / 8, 1)
                                
                                lamp_sets.append({
                                    'F3': F3_seq, 
                                    'FIP': FIP_seq, 
                                    'LoopF': LF_seq,
                                    'LoopB': LB_seq,
                                    'BIP': BIP_seq, 
                                    'B3': B3_seq,
                                    'Start': r1['Location'], 
                                    'End': r6['Location'] + r6['AlignWidth'],
                                    'Ampsize': ampsize, 
                                    'F3_Tm': r1['Tm'], 
                                    'B3_Tm': r6['Tm'],
                                    'Avg_Cons': avg_cons
                                })
                                
                                # Early exit if we have enough sets
                                if len(lamp_sets) >= max_results:
                                    return pd.DataFrame(lamp_sets)
                                    
        return pd.DataFrame(lamp_sets)
