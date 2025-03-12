import pandas as pd
from Bio.Seq import Seq

df=pd.read_csv("primer.csv")

minampsize = int(input("Enter the min amplicon size: "))
maxampsize = int(input("Enter the max amplicon size: "))
minprimertm = int(input("Enter min primer Tm: "))
#minprobetm = int(input("Enter min probe Tm: "))
prdiff = int(input("Max diff betwn fwd and rev primer lengths: "))

listprimern=[]
for n in df.index:
	if float(df["Tm"][n])>=minprimertm:
		listprimern.append(n)

Forward=[];Reverse=[];Start=[];End=[]
TmF=[];TmR=[];GCF=[];GCR=[]
LenF=[];LenR=[];Ampsize=[];#Probe=[]

for item1 in listprimern:
    for item2 in listprimern:
        primerdf=abs(df["Length"][item1]-df["Length"][item2])
        amplicon=(df["Location"][item2]+df["Length"][item2])-(df["Location"][item1])       

        if amplicon >= minampsize and amplicon <= maxampsize and primerdf<=prdiff:
            #for item3 in listproben:
                #if (int(df["PrimerLoc"][item3])>(int(df["PrimerLoc"][item1])+int(df["Length"][item1]))) and (int(df["PrimerLoc"][item3])<(int(df["PrimerLoc"][item2])-int(df["Length"][item3]))):
                    Forward.append(Seq(df["Sequence"][item1]))
                    Reverse.append(Seq(df["Sequence"][item2]).reverse_complement())
                    #Probe.append(Seq(df["PrimerSeq"][item3]))
                    a=df["Location"][item1]
                    b=df["Location"][item2]+df["Length"][item2]-1
                    Start.append(a)
                    End.append(b)
                    TmF.append(float(df["Tm"][item1]))
                    TmR.append(float(df["Tm"][item2]))
                    GCF.append(df["GC"][item1])
                    GCR.append(df["GC"][item2])
                    LenF.append(df["Length"][item1])
                    LenR.append(df["Length"][item2])
                    Ampsize.append(amplicon)
                    
primerpair=pd.DataFrame(list(zip(Forward,Reverse,Start,End,TmF,TmR,GCF,GCR,LenF,LenR,Ampsize)),columns=['Forward','Reverse','Start','End','TmF','TmR','GCF','GCR','LenF','LenR','Ampsize'])

primerpair.to_csv("primerpair"+".csv",index=False)


