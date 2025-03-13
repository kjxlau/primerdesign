<h1>PRIMER DESIGN TOOL</h1>

Hello. I have designed a series of 4-part python scripts to generate candidate primers and probes for Polymerase Chain Reaction (PCR) in your experiments.

I will explain each of the four part python script as follows:

<br />

**Step1_DownloadFasta.py:**

Users will be prompted to enter 

(1) No of fasta entries to download

(2) Organism of interest or NCBI taxa ID

https://www.ncbi.nlm.nih.gov/nuccore/advanced 

Use NCBI advanced search tool to get search strings to download highly specific sequences of interest

Example 1:

>(Human Papillomavirus[Organism]) AND E6[Gene Name] 

Example 2:

>(Dengue Virus[Organism]) AND NS1[Gene Name] 

(3) Filename to save in .fasta format. Users will need to then align the file using either clustalW or mafft multiple sequence aligment programs to generate an aligned sequence fasta file.

Replace the fasta header of the file with the target strain name, so that Step2_FindConsensusSeq.py script can identify which sequences to look at to shortlist primer sequences.

Example 1:

> \>HPV16 (fasta header)

Example 2:

> \>DENV1 (fasta header)

<br />

**Step2_FindConsensusSeq.py:**

Users will be prompted to enter 

(1) Alignment file from step 1

(2) The targeted serotype/species/strain

Example 1:

>HPV16

Example 2:

>DENV1

(3) No of degenerate bases 

(4) The percent cutoff for no of mutations allowed acrossed sequences (recommended: <0.2)

The script will shortlist all potential sequences that are good primer candidates.

<br />

**Step3_PrimerScreening.py:**

Users will be prompted to enter 

(1) Minimum and maximum amplicon size

(2) Melting Temperature of amplicon product

(3) Max difference in length of flanking primer pair

To screen for potential primer sets

<br />

**Step4_ProbeSelect.py:**

Users will be prompted to enter 

(1) Minimum probe melting temperature required

To screen for potential Taqman probe.

It is recommended to set the melting temperature to be at least 5 deg C higher relative to the primer pair.

Please also check for secondary structure (hairpin) using other tools like OligoAnalyzer

https://www.idtdna.com/pages/tools/oligoanalyzer

<br />

Have fun developing your own PCR test kit.

For any queries, you can reach out to me via email at kennyjxlau@gmail.com


<br />
Best regards,

Kenny
