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