import pandas as pd
from Bio.Seq import Seq

# 1. Load data
df = pd.read_csv("primer.csv")

minprimertm = float(input("Enter min primer Tm: "))

# Standard LAMP distance constraints (gap sizes between regions in bp)
# You can tweak these based on your specific design needs!
d_F3_F2_min, d_F3_F2_max = 0, 60
d_F2_F1_min, d_F2_F1_max = 40, 60
d_F1_B1_min, d_F1_B1_max = 10, 60  # Loop region
d_B1_B2_min, d_B1_B2_max = 40, 60
d_B2_B3_min, d_B2_B3_max = 0, 60

print("Filtering and searching for LAMP sets...")

# 2. Filter by Tm and sort by Location 
df_filtered = df[df["Tm"] >= minprimertm].sort_values("Location").reset_index(drop=True)
records = df_filtered.to_dict('records')
n = len(records)

lamp_sets = []

# 3. Search for the 6 consecutive regions
for i in range(n):
    r1 = records[i] # F3 Region
    for j in range(i+1, n):
        r2 = records[j] # F2 Region
        dist1 = r2['Location'] - (r1['Location'] + r1['Length'])
        if dist1 > d_F3_F2_max: break       # Stop searching forward if gap is too big
        if dist1 < d_F3_F2_min: continue    # Skip if gap is too small
        
        for k in range(j+1, n):
            r3 = records[k] # F1 Region
            dist2 = r3['Location'] - (r2['Location'] + r2['Length'])
            if dist2 > d_F2_F1_max: break
            if dist2 < d_F2_F1_min: continue
            
            for l in range(k+1, n):
                r4 = records[l] # B1c Region
                dist3 = r4['Location'] - (r3['Location'] + r3['Length'])
                if dist3 > d_F1_B1_max: break
                if dist3 < d_F1_B1_min: continue
                
                for m in range(l+1, n):
                    r5 = records[m] # B2c Region
                    dist4 = r5['Location'] - (r4['Location'] + r4['Length'])
                    if dist4 > d_B1_B2_max: break
                    if dist4 < d_B1_B2_min: continue
                    
                    for p in range(m+1, n):
                        r6 = records[p] # B3c Region
                        dist5 = r6['Location'] - (r5['Location'] + r5['Length'])
                        if dist5 > d_B2_B3_max: break
                        if dist5 < d_B2_B3_min: continue
                        
                        # -- If we reach here, we found a valid 6-region LAMP set! --
                        
                        # Construct the 4 standard LAMP Primers
                        F3_seq = r1['Sequence']
                        
                        # FIP = F1c + F2 (F1c is reverse complement of F1)
                        F1c = str(Seq(r3['Sequence']).reverse_complement())
                        FIP_seq = F1c + r2['Sequence']
                        
                        # BIP = B1c + B2 (B2 is reverse complement of B2c)
                        B2 = str(Seq(r5['Sequence']).reverse_complement())
                        BIP_seq = r4['Sequence'] + B2
                        
                        # B3 = Reverse complement of B3c
                        B3_seq = str(Seq(r6['Sequence']).reverse_complement())
                        
                        # Calculate total span (Amplicon size from F3 start to B3 end)
                        ampsize = (r6['Location'] + r6['Length']) - r1['Location']
                        
                        lamp_sets.append({
                            'F3': F3_seq,
                            'FIP': FIP_seq,
                            'BIP': BIP_seq,
                            'B3': B3_seq,
                            'Total_Amplicon_Size': ampsize,
                            'Start_Loc': r1['Location'],
                            'End_Loc': r6['Location'] + r6['Length'],
                            'F3_Tm': r1['Tm'],
                            'B3_Tm': r6['Tm']
                        })

# 4. Save to CSV
if lamp_sets:
    primerpair = pd.DataFrame(lamp_sets)
    primerpair.to_csv("lamp_primers.csv", index=False)
    print(f"Success! Found {len(lamp_sets)} potential LAMP primer sets. Saved to 'lamp_primers.csv'")
else:
    print("No valid LAMP primer sets found with the current constraints.")
