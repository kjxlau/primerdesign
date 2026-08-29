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