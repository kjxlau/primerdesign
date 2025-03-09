A series of python scripts to generate primers for Polymerase Chain Reaction.

Step1_DownloadFasta.py: Users will be prompted to enter (1) No of fasta entries to download, (2) Organism of interest or NCBI taxa ID, (3) Countries of interest and (4) Filename to save in .fasta format. Users will need to then align the file using either clustalW or mafft multiple sequence aligment programs to generate an aligned sequence fasta file.

Step2_FindConsensusSeq.py: Users will be prompted to enter (1) Alignment file from step 1, (2) The targeted serotype/species/strain, (3) No of degenerate bases and (4) The percent cutoff for no of mutations allowed acrossed sequences 

Step3_PrimerScreening.py: To screen for potential primer sets

Step4_ProbeSelect.py: To screen for potential Taqman probe.

Kenny
