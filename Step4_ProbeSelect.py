import pandas as pd
from Bio.Seq import Seq

df=pd.read_csv("primer.csv")
df1=pd.read_csv("primerpair.csv")

minprobetm = int(input("Enter min probe Tm: "))

listproben=[]
for n in df.index:
	if float(df["PrimerTm"][n])>=minprobetm:
		listproben.append(n)
		
listprimern=[]
for n in df1.index:
	listprimern.append(n)

Forward=[];Reverse=[];Start=[];End=[]
TmF=[];TmR=[];GCF=[];GCR=[]
LenF=[];LenR=[];Ampsize=[]
Probe=[];PLoc=[];PTm=[]

count=0

for item1 in listproben:
    for item2 in listprimern:
        if (int(df["PrimerLoc"][item1])>(int(df1["Start"][item2])+int(df1["LenF"][item2]))) and (int(df["PrimerLoc"][item1])<(int(df1["End"][item2])-int(df["Length"][item1]))):
            count+=1
            if count<=10:
                Forward.append(Seq(df1["Forward"][item2]))
                Reverse.append(Seq(df1["Reverse"][item2]))
                Probe.append(Seq(df["PrimerSeq"][item1]))
                Start.append(df1["Start"][item2])
                End.append(df1["End"][item2])
                TmF.append(float(df1["TmF"][item2]))
                TmR.append(float(df1["TmR"][item2]))
                GCF.append(df1["GCF"][item2])
                GCR.append(df1["GCR"][item2])
                LenF.append(df1["LenF"][item2])
                LenR.append(df1["LenR"][item2])
                Ampsize.append(df1["Ampsize"][item2])
                PLoc.append(df["PrimerLoc"][item1])
                PTm.append(df["PrimerTm"][item1])
            else:
                break

probeselect=pd.DataFrame(list(zip(Forward,Reverse,Start,End,TmF,TmR,GCF,GCR,LenF,LenR,Ampsize,Probe,PLoc,PTm)),columns=['Forward','Reverse','Start','End','TmF','TmR','GCF','GCR','LenF','LenR','Ampsize','Probe','PLoc','PTm'])

probeselect.to_csv("probeselect"+".csv",index=False)



