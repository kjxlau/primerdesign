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
