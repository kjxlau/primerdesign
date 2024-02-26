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
from pandas import ExcelWriter
from pandas import ExcelFile

#import the functions written in function.py
from Functions import IUB_to_regexp
from Functions import FindSeq

#Step 1: Read the multiple sequence alignment result in
alignment_file=input("The alignment file used to design primers>>>")
alignment=AlignIO.read(alignment_file,'fasta')

#Ask user for target species name or serotype or genotype
target=input("Enter the target species/serotype/genotype: ")

#Step 2: Define the sliding window size to search consensus sequence
minconsensus,maxconsensus=input("Set min and max consensus size>>>").strip().split()
print("Minimum size range from "+minconsensus+"bp to "+maxconsensus+"bp",'\n')
minconsensus=int(minconsensus);maxconsensus=int(maxconsensus)

#Step 3: Set no of degenerate nucleotides allowed in consensus sequence
degenConsensus=int(input('Set the number of degenerate nucleotide for consensus seq: '))
print("Consensus sequence has "+str(degenConsensus)+" degeneracy",'\n')

