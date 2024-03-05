#import biopython package
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from Bio.Seq import Seq
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment as msa
from Bio.SeqUtils import MeltingTemp as mt

#import essential packages
import re
import numpy as np
import math,string,os,sys,shutil,csv,openpyxl,glob
from statistics import mean

#import pandas for dataframe handling
import pandas as pd

#import the functions written in function.py
from Functions import IUB_to_regexp
from Functions import FindSeq

#Step 1: Read the multiple sequence alignment result in
alignment_file=input("The alignment file for primer design: ")
alignment=AlignIO.read(alignment_file,'fasta')

#Ask user for target species name or serotype or genotype
target=input("Enter the target species/serotype/genotype: ")

#Step 2: Define the sliding window size to search consensus sequence
minconsensus=20;maxconsensus=30;

#Step 3: Set no of degenerate nucleotides allowed in consensus sequence
degenConsensus=int(input('Set the number of degenerate nucleotide for consensus seq: '))

#Step 4: Set % cutoff for no of mutations allowed across consensus sequence
Cutoff=float(input("Set percent cutoff of mutant ratio (0 to 1): "))
print("Cutoff ratio is set as "+str(Cutoff),'\n')

#Step 5: Add aligned sequences into array prior to sliding window process
species=[];seqdes=[]
species_align=msa([])
for record in alignment:
	seqdes=re.split('\ |\[|\]|\,|\(|\)|\_|\-',record.description)
	align_array=np.array([list(record.seq) for record in alignment],dtype=str)
	if target in seqdes[0:2]:
        	species.append(record.seq)
	species_align.add_sequence(record.id,str(record.seq))
	spec_array=np.array(species,dtype=str)

summary_align=AlignInfo.SummaryInfo(species_align)
spec_consensus=summary_align.dumb_consensus()

SpecArrayLength=len(spec_array)
print("No of sequences that contains target species/serotype: "+str(SpecArrayLength))
FullLength=len(spec_array.T)
print("Length of aligned sequences: "+str(FullLength))
FullArrayLength=len(align_array)
print("No of entered sequences in alignment file: "+str(FullArrayLength))

#Step 6: Perform the sliding window to search for consensus sequence
seqlist,seqloc,seqGC,seqTm=FindSeq(spec_array,align_array,minconsensus,
	maxconsensus,SpecArrayLength,FullLength,FullArrayLength,degenConsensus,Cutoff)
output=pd.DataFrame(list(zip(seqlist,seqloc,seqGC,seqTm)),
	columns=['Sequence','Location','GC','Tm'])

#Step 7: Remove dashes and find the length of primer sequences
NewPrimerSeq=[];NewPrimerLen=[]
for item in output["Sequence"]:
	seq=str(item).replace('-','')
	length=len(seq)
	NewPrimerSeq.append(seq)
	NewPrimerLen.append(length)

#Step 8: Save the output in a csv file
output=pd.DataFrame(list(zip(NewPrimerSeq,seqloc,NewPrimerLen,seqGC,seqTm)),
	columns=['Sequence','Location','Length','GC','Tm'])
output.to_csv("Lamp_Primers.csv",index=False)
