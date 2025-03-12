Hello. I have designed a series of 4-part python scripts to generate candidate primers and probes for Polymerase Chain Reaction (PCR) in your experiments.

I will explain each of the four part python script as follows:
**Step1_DownloadFasta.py: **
Users will be prompted to enter (1) No of fasta entries to download, (2) Organism of interest or NCBI taxa ID, and (3) Filename to save in .fasta format. Users will need to then align the file using either clustalW or mafft multiple sequence aligment programs to generate an aligned sequence fasta file.

**Step2_FindConsensusSeq.py: **
Users will be prompted to enter (1) Alignment file from step 1, (2) The targeted serotype/species/strain, (3) No of degenerate bases and (4) The percent cutoff for no of mutations allowed acrossed sequences 

**Step3_PrimerScreening.py: **
To screen for potential primer sets

**Step4_ProbeSelect.py: **
To screen for potential Taqman probe.

Have fun developing your own test PCR probes!

For any queries, you can drop me an email at kennyjxlau@gmail.com

Best regards,
Kenny
