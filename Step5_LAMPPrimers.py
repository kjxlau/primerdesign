import pandas as pd
from Bio.Seq import Seq

# 1. Load data
df = pd.read_csv("primer.csv")
minprimertm = float(input("Enter min primer Tm: "))

# Distance constraints (gap sizes between regions in bp)
d_F3_F2_min, d_F3_F2_max = 0, 60
d_F2_F1_min, d_F2_F1_max = 40, 60
d_F1_B1_min, d_F1_B1_max = 10, 60
d_B1_B2_min, d_B1_B2_max = 40, 60
d_B2_B3_min, d_B2_B3_max = 0, 60

print("Filtering and searching for 6-primer LAMP sets (including LF & LB)...")

# 2. Filter by Tm and sort by Location
df_filtered = df[df["Tm"] >= minprimertm].sort_values("Location").reset_index(drop=True)
records = df_filtered.to_dict('records')
n = len(records)
lamp_sets = []

# Helper function to find a loop primer candidate within a region
def find_loop_primer(records, start_bound, end_bound, is_reverse=False):
    for r in records:
        r_start = r['Location']
        r_end = r['Location'] + r['Length']
        
        # Check if the candidate fits entirely within the loop boundaries
        if r_start >= start_bound and r_end <= end_bound:
            seq = str(Seq(r['Sequence']).reverse_complement()) if is_reverse else r['Sequence']
            return seq, r['Tm']
    return None, None

# 3. Search for consecutive regions
for i in range(n):
    r1 = records[i]  # F3
    
    for j in range(i+1, n):
        r2 = records[j]  # F2
        dist1 = r2['Location'] - (r1['Location'] + r1['Length'])
        if dist1 > d_F3_F2_max: break
        if dist1 < d_F3_F2_min: continue
        
        for k in range(j+1, n):
            r3 = records[k]  # F1
            dist2 = r3['Location'] - (r2['Location'] + r2['Length'])
            if dist2 > d_F2_F1_max: break
            if dist2 < d_F2_F1_min: continue
            
            for l in range(k+1, n):
                r4 = records[l]  # B1c
                dist3 = r4['Location'] - (r3['Location'] + r3['Length'])
                if dist3 > d_F1_B1_max: break
                if dist3 < d_F1_B1_min: continue
                
                for m in range(l+1, n):
                    r5 = records[m]  # B2c
                    dist4 = r5['Location'] - (r4['Location'] + r4['Length'])
                    if dist4 > d_B1_B2_max: break
                    if dist4 < d_B1_B2_min: continue
                    
                    for p in range(m+1, n):
                        r6 = records[p]  # B3c
                        dist5 = r6['Location'] - (r5['Location'] + r5['Length'])
                        if dist5 > d_B2_B3_max: break
                        if dist5 < d_B2_B3_min: continue

                        # 4. Search for LF (between F2 end and F1 start)
                        lf_start_bound = r2['Location'] + r2['Length']
                        lf_end_bound = r3['Location']
                        LF_seq, LF_tm = find_loop_primer(records, lf_start_bound, lf_end_bound, is_reverse=True)

                        # 5. Search for LB (between B1c end and B2c start)
                        lb_start_bound = r4['Location'] + r4['Length']
                        lb_end_bound = r5['Location']
                        LB_seq, LB_tm = find_loop_primer(records, lb_start_bound, lb_end_bound, is_reverse=False)

                        # Construct core primers
                        F3_seq = r1['Sequence']
                        F1c = str(Seq(r3['Sequence']).reverse_complement())
                        FIP_seq = F1c + r2['Sequence']
                        
                        B2 = str(Seq(r5['Sequence']).reverse_complement())
                        BIP_seq = r4['Sequence'] + B2
                        B3_seq = str(Seq(r6['Sequence']).reverse_complement())

                        ampsize = (r6['Location'] + r6['Length']) - r1['Location']

                        lamp_sets.append({
                            'F3': F3_seq,
                            'FIP': FIP_seq,
                            'BIP': BIP_seq,
                            'B3': B3_seq,
                            'LF': LF_seq if LF_seq else "N/A",
                            'LB': LB_seq if LB_seq else "N/A",
                            'LF_Tm': LF_tm if LF_tm else "N/A",
                            'LB_Tm': LB_tm if LB_tm else "N/A",
                            'Total_Amplicon_Size': ampsize,
                            'Start_Loc': r1['Location'],
                            'End_Loc': r6['Location'] + r6['Length'],
                            'F3_Tm': r1['Tm'],
                            'B3_Tm': r6['Tm']
                        })

# Save to CSV
if lamp_sets:
    primerpair = pd.DataFrame(lamp_sets)
    primerpair.to_csv("lamp_primers.csv", index=False)
    print(f"Success! Found {len(lamp_sets)} potential LAMP primer sets. Saved to 'lamp_primers.csv'")
else:
    print("No valid LAMP primer sets found with the current constraints.")
