<h1>PRIMER DESIGN TOOL</h1>

Hello. I have designed a series of 4-part python scripts to generate candidate primers and probes for Polymerase Chain Reaction (PCR) in your experiments.

I will explain each of the four part python script as follows:

<br />

**Step1_DownloadFasta.py:**

Users will be prompted to enter 

(1) No of fasta entries to download

(2) Organism of interest or NCBI taxa ID

eg. Human Papillomavirus [Organism] and E6 [Gene]

(3) Filename to save in .fasta format. Users will need to then align the file using either clustalW or mafft multiple sequence aligment programs to generate an aligned sequence fasta file.

Label the fasta header of the file with the target strain name, so that Step2_FindConsensusSeq.py can identify which sequences to look at to shortlist primer sequences.

<br />

**Step2_FindConsensusSeq.py:**

Users will be prompted to enter 

(1) Alignment file from step 1

(2) The targeted serotype/species/strain

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

Have fun developing your own PCR test kit!

For any queries, you can drop me an email at kennyjxlau@gmail.com


<br />
Best regards,

Kenny
