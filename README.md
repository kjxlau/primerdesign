# 🧬 Automated Multi-Agent Bioinformatics Pipeline

An AI-orchestrated bioinformatics workflow for automated **qPCR** and **LAMP** primer and probe design. 

Built using **LangGraph**, **LangChain**, and **Biopython**, this tool automates the tedious and error-prone process of manual primer design. By stringing together specialized "agents," the pipeline autonomously fetches sequences, aligns them, identifies conserved regions, calculates thermodynamic properties, and outputs highly specific primer/probe sets for infectious disease diagnostics or genetic research.

---

## 🔬 1. Overview

### What are we trying to solve?
Designing robust qPCR and LAMP assays traditionally requires hopping between multiple databases (NCBI), alignment tools (Clustal/MAFFT), and primer design software (Primer3), followed by manual thermodynamic calculations. This manual process is time-consuming and prone to human error.

### The Solution
This pipeline orchestrates a **multi-agent state graph** to automate the workflow from end to end. You simply input a target organism and gene, and the system handles the rest.

### ⚙️ Pipeline Workflow
The workflow runs as a directed graph with the following sequential and parallel nodes:
1. **🔍 Search Agent (`search`):** Connects to the NCBI Nucleotide database via Entrez to smartly query and fetch relevant FASTA sequences (e.g., 16S, ITS, or specific genes) while filtering out incomplete records.
2. **📏 Alignment Agent (`align`):** Processes raw sequences and performs multiple sequence alignment (MSA).
3. **📊 Analyst Agent (`analyze`):** Analyzes the alignment to identify highly conserved regions suitable for primer targeting.
4. **🧪 qPCR Agent (`qpcr`):** *(Runs in parallel)* Screens for forward/reverse primer pairs and selects an internal fluorescent probe based on strict thermodynamic parameters (Tm, amplicon size, GC content).
5. **💡 LAMP Agent (`lamp`):** *(Runs in parallel)* Designs Loop-Mediated Isothermal Amplification primer sets (F3, B3, FIP, BIP, LoopF, LoopB) for rapid diagnostics.
6. **🤖 Reporter Agent (`report`):** Uses **OpenAI (GPT-4o)** to analyze the run context and generate a final human-readable summary report.

---

## 🛠️ 2. Setup Instructions

### Prerequisites
* **Python 3.8+**
* An **NCBI Account** (for API key and email)
* An **OpenAI API Key** (for the reporting agent)

### Step 1: Clone the Repository
```bash
git clone https://github.com/yourusername/bioinformatics-pipeline.git
cd bioinformatics-pipeline
```

### Step 2: Create a Virtual Environment (Recommended)
```bash
python -m venv venv

# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate
```

### Step 3: Install Dependencies
Ensure you install the required Python libraries:
```bash
pip install biopython pandas langgraph langchain-openai langchain-core python-dotenv requests
```
*(If your AlignmentAgent relies on local tools like MAFFT or ClustalW, ensure they are installed on your system's PATH).*

### Step 4: Configure Environment Variables
In the root directory of the project, create a file named `.env`. 
Add your API keys and credentials to this file:

```env
# .env file
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=your_ncbi_api_key_here
OPENAI_API_KEY=sk-your_openai_api_key_here
```
> **Note:** Never commit your `.env` file to version control. Ensure it is included in your `.gitignore`.

---

## 🚀 3. How to Run the Pipeline

Execute the main script from your terminal. The script is interactive and will prompt you for your targets.

```bash
python main.py
```

### Example Interactive Run:
```text
--- Primer & Probe Design Pipeline Setup ---
Enter organism of interest [default: Escherichia coli]: Mycobacterium tuberculosis
Enter target gene/region [default: 16S]: rpoB
Enter matching keyword [default: Mycobacterium]: Mycobacterium
Enter number of sequences to fetch [default: 10]: 15

🚀 Running pipeline for: 'Mycobacterium tuberculosis' (Gene: 'rpoB')...
```

### 📂 Output Files
Once the pipeline finishes executing, it will output the following in your working directory:
1. **`{organism}_{gene}_alignment.fasta`** - The aligned sequences used for analysis.
2. **`{organism}_{gene}_probe_sets.csv`** - The final filtered list of viable qPCR primers and probes.
3. **`{organism}_{gene}_lamp_sets.csv`** - The final filtered list of LAMP primer sets.
4. **Console Output** - An LLM-generated executive summary detailing the success and results of the pipeline run.