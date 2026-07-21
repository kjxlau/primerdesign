import os
import ollama
import json
from Multiagent_orchestrator import MasterOrchestrator # Importing your previous logic

# --- THE OLLAMA BRAIN ---
class OllamaDirector:
    def __init__(self, model="llama3"):
        self.model = model

    def interpret_user_request(self, user_input):
        """
        Uses Ollama to turn a sentence into a JSON of biological parameters.
        """
        prompt = f"""
        You are a Senior Bioinformatician. Extract parameters from this request: "{user_input}"
        Return ONLY a JSON object with these keys: 
        organism, gene, count, target_kw, min_amp, max_amp, min_tm, min_probe_tm.
        If a value is missing, use these defaults: 
        count: 10, target_kw: same as organism, min_amp: 100, max_amp: 300, min_tm: 55, min_probe_tm: 65.
        
        Example Output: {{"organism": "Bacillus subtilis", "gene": "16S", "count": 10, ...}}
        """
        
        print(f"[OllamaDirector]: Thinking...")
        response = ollama.generate(model=self.model, prompt=prompt)
        
        # Clean the response string to ensure it's valid JSON
        json_str = response['response'].strip()
        # Find the first '{' and last '}' in case the model added conversational text
        json_str = json_str[json_str.find("{"):json_str.rfind("}")+1]
        
        return json.loads(json_str)

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
            'min_amp': params.get('min_amp', 100),
            'max_amp': params.get('max_amp', 300),
            'min_tm': params.get('min_tm', 55),
            'min_probe_tm': params.get('min_probe_tm', 60),
            'max_diff': 5
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
