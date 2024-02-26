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
from SeqCounting_withdegeneracy import IUB_to_regexp
from Functions import FindSeq

#Step1: Read the multiple sequence alignment result in
alignment_file=input("The alignment file used to design primers>>>")
alignment=AlignIO.read(alignment_file,'fasta')
