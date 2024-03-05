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

F3_B3=pd.read_csv("F3_B3.csv")
F2_B2=pd.read_csv("F2_B2.csv")
F1c_B1c=pd.read_csv("F1c_B1c.csv")
LF_LB=pd.read_csv("LF_LB.csv")

listn=[]
for n in LF_LB.index:
	listn.append(n)

F3=[];F3Loc=[];F3Len=[]
B3=[];B3Loc=[];B3Len=[]
FIP=[];FIPLoc=[];FIPLen=[]
BIP=[];BIPLoc=[];BIPLen=[]
LF=[];LFLoc=[];LFLen=[]
LR=[];LRLoc=[];LRLen=[]
for n1 in listn:
    for n2 in listn:
        if int(F3_B3["F3Loc"][n1]+30)<F2_B2["F2Loc"][n2] and
        LF_LB["LFLoc"][n2]<F2_B2["F2Loc"][n2] and
        F2_B2["F2Loc"][n2]<F1c_Bc1["F1cLoc"][n2] and
        F3_B3["B3Loc"][n2]>int(F2_B2["B2Loc"][n1]+30) and
        LF_LB["LBLoc"][n2]<F2_B2["B2Loc"][n2] and
        F2_B2["B2Loc"][n2]<F1c_B1c["B1cLoc"][n2]:
            F3.append(F3_B3["F3"][n1])
            F3Loc.append(F3_B3["F3Loc"][n1])
            F3Len.append(F3_B3["F3Len"][n1])
            FIP.append((F1c_B1c["F1c"][n2])+(LF_LB["LF"][n2])+(F2_B2["F2"][n2]))
            FIPLoc.append(F2_B2["F2Loc"][n2])
            FIPLen.append(F2_B2["F2Len"][n2])
            LF.append(LF_LB["LF"][n2])
            LFLoc.append(LF_LB["LFLoc"][n2])
            LFLen.append(LF_LB["LFLen"][n2])
            B3.append(F3_B3["B3"][n2])
            B3Loc.append(F3_B3["B3Loc"][n2])
            B3Len.append(F3_B3["B3Len"][n2])
            BIP.append((F1c_Bc1["B1c"][n1])+(LF_LB["LB"][n1])+(F2_B2["B2"][n1]))
            BIPLoc.append(F2_B2["B2Loc"][n1])
            BIPLen.append(F2_B2["B2Len"][n1])
            LB.append(LF_LB["LB"][n1])
            LBLoc.append(LF_LB["LBLoc"][n1])
            LBLen.append(LF_LB["LBLen"][n1])
            output=pd.DataFrame(list(zip(F3,F3Loc,F3Len,FIP,FIPLoc,FIPLen,LF,LFLoc,LFLen,LB,LBLoc,LBLen,BIP,BIPLoc,BIPLen,B3,B3Loc,B3Len)),columns=['F3','F3Loc','F3Len','FIP','FIPLoc','FIPLen','LF','LFLoc','LFLen','LB','LBLoc','LBLen','BIP','BIPLoc','BIPLen','B3','B3Loc','B3Len'])
            output.to_csv("final_lamp_primers"+".csv",index=False)
