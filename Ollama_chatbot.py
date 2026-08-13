import os
import json
import ollama
import pandas as pd
from Multiagent_orchestrator import MasterOrchestrator 

# --- THE OLLAMA BRAIN ---
class OllamaDirector:
    def __init__(self, model="llama3"):
        self.model = model

    def interpret_user_request(self, user_input):
        """
        Uses Ollama to turn a natural language sentence into a JSON of biological parameters.
        Now understands requests for both TaqMan and LAMP assays.
        """
        prompt = f"""
        You are a strict biological data extraction system. 
        Extract parameters from this request: "{user_input}"
        Return ONLY a valid JSON object. 
        Do NOT include any markdown formatting like ```json.
        Do NOT include any conversational text.
        All keys and string values MUST be enclosed in strict DOUBLE QUOTES. No trailing commas.
        
        Keys to include: 
        "organism", "gene", "count", "target_kw", "min_amp", "max_amp", "min_tm", "min_probe_tm".
        
        Rules: 
        - If 'count' is missing, default to 15. It MUST NEVER be less than 5.
        - 'target_kw' defaults to the organism name (first word, lower case, e.g., 'bacillus').
        - Default PCR params: min_amp: 100, max_amp: 300, min_tm: 55, min_probe_tm: 60.
        - LAMP relies on min_tm, so ensure min_tm is always provided.
        
        Example Output: 
        {{"organism": "Bacillus subtilis", "gene": "16S rRNA", "count": 15, "target_kw": "subtilis", "min_amp": 100, "max_amp": 300, "min_tm": 55, "min_probe_tm": 60}}
        """
        
        print(f"[OllamaDirector]: Thinking and structuring parameters...")
        response = ollama.generate(model=self.model, prompt=prompt)
        
        json_str = response['response'].strip()
        
        try:
            # Safely extract block if LLM adds text around the JSON
            if "{" in json_str and "}" in json_str:
                json_str = json_str[json_str.find("{"):json_str.rfind("}")+1]
            
            data = json.loads(json_str)
            
            # Guarantee all keys exist
            return {
                "organism": data.get("organism", "Unknown_Org"),
                "gene": data.get("gene", "Unknown_Gene"),
                "count": data.get("count", 15),
                "target_kw": data.get("target_kw", data.get("organism", "").split()[0].lower()),
                "min_amp": data.get("min_amp", 100),
                "max_amp": data.get("max_amp", 300),
                "min_tm": data.get("min_tm", 55),
                "min_probe_tm": data.get("min_probe_tm", 60)
            }
        except Exception as e:
            print(f"[OllamaDirector]: Error parsing JSON - {e}")
            print(f"[OllamaDirector]: Raw LLM output was:\n{response['response']}")
            return None
            
    def explain_results(self, qpcr_file, lamp_file):
        """
        Ollama reads the final results (qPCR and/or LAMP) and gives you a summary.
        """
        prompt = "You are a senior molecular biologist. Here are the top assays discovered by the automated pipeline:\n\n"
        
        has_results = False
        
        # 1. Inject qPCR Results if they exist
        if os.path.exists(qpcr_file):
            df_q = pd.read_csv(qpcr_file).head(2)
            prompt += f"--- Standard TaqMan qPCR Assays ---\n{df_q.to_string()}\n\n"
            has_results = True
            
        # 2. Inject LAMP Results if they exist
        if os.path.exists(lamp_file):
            df_l = pd.read_csv(lamp_file).head(2)
            prompt += f"--- Isothermal LAMP Assays ---\n{df_l.to_string()}\n\n"
            has_results = True
            
        if not has_results:
            return "No results files were generated. The pipeline may have failed to find highly conserved regions."

        prompt += """
        Summarize these results in a concise, professional manner. 
        Explain why these properties (Tm, Amplicon Size, Conservation, or spatial layout for LAMP) make them good candidates for experimental validation.
        """
        
        print(f"[OllamaDirector]: Analyzing the generated assay data...")
        response = ollama.generate(model=self.model, prompt=prompt)
        return response['response']

# --- INTEGRATED EXECUTION ---
if __name__ == "__main__":
    # 1. Initialize the LLM Director
    director = OllamaDirector(model="llama3")
    
    # 2. Get Natural Language input
    print("\n" + "="*50)
    print(" 🧬 OLLAMA MULTI-AGENT BIO-SYSTEM (TaqMan & LAMP) ")
    print("="*50)
    user_req = input("What do you want to design today?\n(e.g., 'Design both a TaqMan and LAMP assay for ITS in Candida albicans')\n> ")
    
    # 3. Ollama extracts the logic
    try:
        params = director.interpret_user_request(user_req)
        if not params:
            exit(1)
            
        print(f"\n[Ollama Extracted Params]: {params}\n")
        
        # 4. Trigger the Python Orchestrator
        pcr_params = {
            'min_amp': params.get('min_amp', 100),
            'max_amp': params.get('max_amp', 300),
            'min_tm': params.get('min_tm', 55),
            'min_probe_tm': params.get('min_probe_tm', 60),
            'max_diff': 2
        }
        
        # Instantiate and run the Orchestrator
        bio_master = MasterOrchestrator()
        bio_master.execute(
            params['organism'], 
            params['gene'], 
            params['count'], 
            params['target_kw'], 
            pcr_params
        )
        
        # 5. Determine filenames based on Orchestrator logic
        base_name = f"{params['organism']}_{params['gene']}".replace(" ","_").lower()
        qpcr_filename = f"{base_name}_probe_sets.csv"
        lamp_filename = f"{base_name}_lamp_sets.csv"
        
        # 6. Ollama explains the output
        print("\n" + "="*50)
        summary = director.explain_results(qpcr_filename, lamp_filename)
        print(" OLLAMA POST-ANALYSIS ")
        print("="*50)
        print(summary)
            
    except Exception as e:
        print(f"Workflow Error: {e}")
