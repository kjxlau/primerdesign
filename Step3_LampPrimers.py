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

#Step 1: Import primer sequence above Tm of 60
df=pd.read_csv("Lamp_Primers1.csv")
listprimern=[]
for n in df.index:
	if float(df["Tm"][n])>=60:
		listprimern.append(n)

F3=[];F3Loc=[];F3Len=[]
B3=[];B3Loc=[];B3Len=[]
count=0

#Step 2: Find F3 and B3 primers
for n1 in listprimern:
    for n2 in listprimern:
        primerdf=abs(df["Location"][n1]-df["Location"][n2])
        if (primerdf<=400) and int(df["Location"][n1]+df["Length"][n1]+250)<int(df["Location"][n2]):
            count+=1
            if count<=1:
                F3.append(Seq(df["Sequence"][n1]))
                F3Loc.append(df["Location"][n1])
                F3Len.append(df["Length"][n1])
                B3.append(Seq(df["Sequence"][n2]).reverse_complement())
                B3Loc.append(df["Location"][n2])
                B3Len.append(df["Length"][n2])
                df2=pd.DataFrame(list(zip(F3,F3Loc,F3Len,B3,B3Loc,B3Len)),columns=['F3','F3Loc','F3Len','B3','B3Loc','B3Len'])
                df2.to_csv("F3_B3"+".csv",index=False)
        else:
            count=0;

F2=[];F2Loc=[];F2Len=[]
B2=[];B2Loc=[];B2Len=[]
count=0
#Find F2 and B2 primers
listprimern2=[]
for n in df2.index:
	listprimern2.append(n)
for n1 in listprimern2:
    for n2 in listprimern2:
        primerdf2=abs(df2["F3Loc"][n1]-df2["B3Loc"][n2])
        if (primerdf2<=160 and
        int(df2["F3Loc"][n2]+df2["F3Len"][n2]+40)<int(df2["F3Loc"][n1]) and
        int(df2["B3Loc"][n1]-40)>int(df2["B3Loc"][n2])):
            count+=1
            if count<=1:
                F2.append(df2["F3"][n1])
                F2Loc.append(df2["F3Loc"][n1])
                F2Len.append(df2["F3Len"][n1])
                B2.append(df2["B3"][n2])
                B2Loc.append(df2["B3Loc"][n2])
                B2Len.append(df2["B3Len"][n2])
                df3=pd.DataFrame(list(zip(F2,F2Loc,F2Len,B2,B2Loc,B2Len)),columns=['F2','F2Loc','F2Len','B2','B2Loc','B2Len'])
                df3.to_csv("F2_B2"+".csv",index=False)
        else:
            count=0;
            
F1c=[];F1cLoc=[];F1cLen=[]
B1c=[];B1cLoc=[];B1cLen=[]
count=0
#Find F1c and B1c primers
listprimern3=[]
for n in df3.index:
	listprimern3.append(n)
for n1 in listprimern3:
    for n2 in listprimern3:
        if int(df3["F2Loc"][n2]+df3["F2Len"][n2]+40)<int(df3["F2Loc"][n1]) and (int(df3["B2Loc"][n2]-40)>int(df3["B2Loc"][n1])) :
            count+=1
            if count<=1:
                F1c.append((df3["F2"][n1]).reverse_complement())
                F1cLoc.append(df3["F2Loc"][n1])
                F1cLen.append(df3["F2Len"][n1])
                B1c.append((df3["B2"][n2]).reverse_complement())
                B1cLoc.append(df2["B2Loc"][n2])
                B1cLen.append(df2["B2Len"][n2])
                df4=pd.DataFrame(list(zip(F1c,F1cLoc,F1cLen,B1c,B1cLoc,B1cLen)),columns=['F1c','F1cLoc','F1cLen','B1c','B1cLoc','B1cLen'])
                df4.to_csv("F1c_B1c"+".csv",index=False)
        else:
            count=0;

