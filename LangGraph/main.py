import os
import time
import requests
from io import StringIO
from typing import TypedDict, Optional
import pandas as pd
from dotenv import load_dotenv
from langgraph.graph import StateGraph, START, END

from Bio import Entrez, AlignIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Align import AlignInfo, MultipleSeqAlignment

from langchain_openai import ChatOpenAI
from langchain_core.messages import SystemMessage, HumanMessage

# --- LOAD ENVIRONMENT VARIABLES ---
load_dotenv() 

# --- IMPORT AGENTS ---
from agents.search_agent import SearchAgent
from agents.alignment_agent import AlignmentAgent
from agents.analyst_agent import AnalystAgent
from agents.screening_agent import ScreeningAgent
from agents.probe_agent import ProbeAgent
from agents.lamp_agent import LampAgent

# --- INSTANTIATE AGENTS ---
searcher = SearchAgent(os.getenv("NCBI_EMAIL"), os.getenv("NCBI_API_KEY"))
aligner = AlignmentAgent(os.getenv("NCBI_EMAIL"))
analyst = AnalystAgent()
screener = ScreeningAgent()
probe_agent = ProbeAgent()
lamp_agent = LampAgent()

# --- STATE SCHEMA ---
class PipelineState(TypedDict):
    organism: str
    gene: str
    count: int
    kw: str
    pcr_params: dict
    raw_sequences: Optional[str]
    alignment_fasta: Optional[str]
    candidates: Optional[pd.DataFrame]
    report: Optional[str]

# --- PIPELINE NODES ---
def search_node(state: PipelineState) -> dict:
    count = state["count"] if isinstance(state["count"], int) and state["count"] >= 2 else 10
    raw = searcher.fetch_sequences(state["organism"], state["gene"], count)
    return {"raw_sequences": raw}

def align_node(state: PipelineState) -> dict:
    if not state.get("raw_sequences") or state["raw_sequences"].count('>') < 2:
        return {"alignment_fasta": None}
    aln = aligner.align(state["raw_sequences"])
    if aln:
        fn = f"{state['organism']}_{state['gene']}_alignment.fasta".replace(" ", "_").lower()
        with open(fn, "w") as f:
            f.write(aln)
    return {"alignment_fasta": aln}

def analyze_node(state: PipelineState) -> dict:
    if not state.get("alignment_fasta"):
        return {"candidates": None}
    candidates = analyst.calculate_stats(state["alignment_fasta"], state["kw"], 2, 0.05)
    return {"candidates": candidates if not candidates.empty else None}

def qpcr_node(state: PipelineState) -> dict:
    candidates = state.get("candidates")
    pcr_params = state["pcr_params"]
    if candidates is not None:
        pairs = screener.screen_pairs(candidates, pcr_params['min_amp'], pcr_params['max_amp'], pcr_params['min_tm'], pcr_params['max_diff'])
        if not pairs.empty:
            final_probe_df = probe_agent.select_probes(candidates, pairs, pcr_params['min_probe_tm'])
            if not final_probe_df.empty:
                fn = f"{state['organism']}_{state['gene']}_probe_sets.csv".replace(" ", "_").lower()
                final_probe_df.to_csv(fn, index=False)
    return {}

def lamp_node(state: PipelineState) -> dict:
    candidates = state.get("candidates")
    pcr_params = state["pcr_params"]
    if candidates is not None:
        lamp_df = lamp_agent.design_lamp(candidates, pcr_params['min_tm'])
        if not lamp_df.empty:
            fn = f"{state['organism']}_{state['gene']}_lamp_sets.csv".replace(" ", "_").lower()
            lamp_df.to_csv(fn, index=False)
    return {}

def report_node(state: PipelineState) -> dict:
    """Uses OpenAI to write a final summary report of the bioinformatics run."""
    llm = ChatOpenAI(model="gpt-4o", temperature=0.2)
    
    organism = state.get("organism")
    gene = state.get("gene")
    candidates_found = state.get("candidates") is not None
    
    # Fix the LLM prompt so it doesn't hallucinate
    csv_status = "CSV files were successfully generated." if candidates_found else "No CSV files were generated because no viable candidates were found."
    
    sys_msg = SystemMessage(content="You are an expert bioinformatics assistant.")
    human_msg = HumanMessage(content=(
        f"A primer design pipeline has just finished for the organism '{organism}' and target gene '{gene}'.\n"
        f"Were viable alignment candidates found? {'Yes' if candidates_found else 'No'}.\n"
        f"Resulting file status: {csv_status}\n"
        f"Please write a brief, professional 3-sentence summary stating that the automated qPCR and "
        f"LAMP primer design process has concluded. Accurately reflect the file generation status provided."
    ))
    
    response = llm.invoke([sys_msg, human_msg])
    
    print("\n" + "="*40)
    print("🤖 OPENAI FINAL REPORT:")
    print("="*40)
    print(response.content)
    print("="*40 + "\n")
    
    return {"report": response.content}

# --- BUILD AND COMPILE STATEGRAPH ---
builder = StateGraph(PipelineState)

# Add Nodes
builder.add_node("search", search_node)
builder.add_node("align", align_node)
builder.add_node("analyze", analyze_node)
builder.add_node("qpcr", qpcr_node)
builder.add_node("lamp", lamp_node)
builder.add_node("report", report_node)

# Add Sequential & Parallel Edges
builder.add_edge(START, "search")
builder.add_edge("search", "align")
builder.add_edge("align", "analyze")
builder.add_edge("analyze", "qpcr")
builder.add_edge("analyze", "lamp")
builder.add_edge(["qpcr", "lamp"], "report")
builder.add_edge("report", END)

app = builder.compile()

# --- EXECUTION ---
if __name__ == "__main__":
    print("\n" + "="*45)
    print("🧬 BIOINFORMATICS PRIMER DESIGN PIPELINE")
    print("="*45)
    
    # Prompt the user for organism and related parameters
    user_organism = input("Enter organism of interest [default: Escherichia coli]: ").strip()
    organism = user_organism if user_organism else "Escherichia coli"
    
    user_gene = input("Enter target gene/region [default: 16S]: ").strip()
    gene = user_gene if user_gene else "16S"
    
    # Auto-suggest the first word as the genus filter keyword
    default_kw = organism.split()[0] if organism else "Escherichia"
    user_kw = input(f"Enter matching keyword [default: {default_kw}]: ").strip()
    kw = user_kw if user_kw else default_kw
    
    user_count = input("Enter number of sequences to fetch [default: 10]: ").strip()
    count = int(user_count) if user_count.isdigit() else 10

    initial_input = {
        "organism": organism,
        "gene": gene,
        "count": count,
        "kw": kw,
        "pcr_params": {
            'min_amp': 70, 
            'max_amp': 200, 
            'min_tm': 58, 
            'min_probe_tm': 68, 
            'max_diff': 2
        }
    }
    
    print(f"\n🚀 Running pipeline for: '{organism}' (Gene: '{gene}')...\n")
    app.invoke(initial_input)
