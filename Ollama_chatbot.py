import os
import ollama
import json
from Multiagent_orchestrator import MasterOrchestrator 

# --- THE OLLAMA BRAIN ---
class OllamaDirector:
    def __init__(self, model="llama3"):
        self.model = model

    def interpret_user_request(self, user_input):
        """
        Uses Ollama to turn a sentence into a JSON of biological parameters.
        """
        prompt = f"""
        You are a strict data extraction system. Extract parameters from this request: "{user_input}"
        Return ONLY a valid JSON object. 
        Do NOT include any markdown formatting like ```json.
        Do NOT include any conversational text.
        All keys and string values MUST be enclosed in strict DOUBLE QUOTES. No trailing commas.
        
        Keys to include: 
        "organism", "gene", "count", "target_kw", "min_amp", "max_amp", "min_tm", "min_probe_tm".
        
        Rules: 
        - If 'count' is missing, default to 15. It MUST NEVER be less than 5.
        - 'target_kw' defaults to the organism name.
        - Default PCR params: min_amp: 100, max_amp: 300, min_tm: 55, min_probe_tm: 60.
        
        Example Output: 
        {{"organism": "Bacillus subtilis", "gene": "16S rRNA", "count": 15, "target_kw": "Bacillus subtilis", "min_amp": 100, "max_amp": 300, "min_tm": 55, "min_probe_tm": 60}}
        """
        
        print(f"[OllamaDirector]: Thinking...")
        response = ollama.generate(model=self.model, prompt=prompt)
        
        json_str = response['response'].strip()
        
        try:
            # Safely extract block
            if "{" in json_str and "}" in json_str:
                json_str = json_str[json_str.find("{"):json_str.rfind("}")+1]
            
            data = json.loads(json_str)
            
            # Guarantee all keys exist
            return {
                "organism": data.get("organism", "Unknown_Org"),
                "gene": data.get("gene", "Unknown_Gene"),
                "count": data.get("count", 15),
                "target_kw": data.get("target_kw", data.get("organism", "")),
                "min_amp": data.get("min_amp", 100),
                "max_amp": data.get("max_amp", 300),
                "min_tm": data.get("min_tm", 55),
                "min_probe_tm": data.get("min_probe_tm", 60)
            }
        except Exception as e:
            print(f"[OllamaDirector]: Error parsing JSON - {e}")
            print(f"[OllamaDirector]: Raw LLM output was:\n{response['response']}")
            return None
            
    def explain_results(self, csv_file):
        """
        Ollama reads the final result and gives you a summary.
        """
        import pandas as pd
        df = pd.read_csv(csv_file).head(3)
        
        prompt = f"""
        Here are the top 3 qPCR primer/probe sets discovered by the pipeline:
        {df.to_string()}
        
        Summarize these results and explain why they are good for a TaqMan assay.
        """
        response = ollama.generate(model=self.model, prompt=prompt)
        return response['response']

# --- INTEGRATED EXECUTION ---
if __name__ == "__main__":
    # 1. Initialize the LLM Director
    director = OllamaDirector(model="llama3")
    
    # 2. Get Natural Language input
    print("\n--- OLLAMA BIO-SYSTEM ---")
    user_req = input("What do you want to design today?\n(e.g., 'Design a TaqMan assay for 16S RNA in Bacillus subtilis')\n> ")
    
    # 3. Ollama extracts the logic
    try:
        params = director.interpret_user_request(user_req)
        print(f"\n[Ollama] Parameters Extracted: {params}")
        
        # 4. Trigger the Python Orchestrator
        pcr_params = {
            'min_amp': params.get('min_amp') or 100,
            'max_amp': params.get('max_amp') or 300,
            'min_tm': params.get('min_tm') or 55,
            'min_probe_tm': params.get('min_probe_tm') or 60,
            'max_diff': 2
        }
        
        bio_master = MasterOrchestrator()
        bio_master.execute(
            params['organism'], 
            params['gene'], 
            params['count'], 
            params['target_kw'], 
            pcr_params
        )
        
        # 5. Ollama explains the output
        filename = f"{params['organism']}_{params['gene']}_probe_sets.csv".replace(" ","_").lower()
        if os.path.exists(filename):
            summary = director.explain_results(filename)
            print("\n--- OLLAMA ANALYSIS ---")
            print(summary)
            
    except Exception as e:
        print(f"Workflow Error: {e}")
