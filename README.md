<h1>PRIMER DESIGN TOOL</h1>

Hello. I have designed a series of 4-part python scripts to generate candidate primers and probes for Polymerase Chain Reaction (PCR) in your experiments.

I will explain each of the four part python script as follows:

<br />

**Step1_DownloadFasta.py:**

Users will be prompted to enter 

(1) No of fasta entries to download

(2) Organism of interest or NCBI taxa ID

(3) Filename to save in .fasta format. Users will need to then align the file using either clustalW or mafft multiple sequence aligment programs to generate an aligned sequence fasta file.

<br />

**Step2_FindConsensusSeq.py:**

Users will be prompted to enter 

(1) Alignment file from step 1

(2) The targeted serotype/species/strain

(3) No of degenerate bases 

(4) The percent cutoff for no of mutations allowed acrossed sequences 

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

<br />

Have fun developing your own test PCR probes!

For any queries, you can drop me an email at kennyjxlau@gmail.com


<br />
Best regards,

Kenny
