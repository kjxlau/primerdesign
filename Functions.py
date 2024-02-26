#import biopython package
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

#import essential packages
import re
import numpy as np
import math,string,os,sys,shutil,csv,openpyxl,glob
from statistics import mean

#import pandas for dataframe handling
import pandas as pd

from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

def IUB_to_regexp(iub):
 
    regular_expression = ''
    iub2character_class = {
     
        '-' : '-',
        'A' : 'A',
        'C' : 'C',
        'G' : 'G',
        'T' : 'T',
        'R' : '[GA]',
        'Y' : '[CT]',
        'M' : '[AC]',
        'K' : '[GT]',
        'S' : '[GC]',
        'W' : '[AT]',
    }

    for i in range(len(iub)):
    	regular_expression+=iub2character_class[iub[i]]
    return regular_expression

def ConsensusSeq(array,ratio):
	consensus=[]
	for nt in range(len(array.T)):
		A=0;T=0;G=0;C=0;D=0;NuType=0;CountType=0
		for row in range(len(array)):
			if array[row][nt]in{'a','A'}:
				A+=1
			elif array[row][nt]in{'t','T'}:
				T+=1
			elif array[row][nt]in{'c','C'}:
				C+=1
			elif array[row][nt]in{'g','G'}:
				G+=1
			else:
				D+=1
		
		if A<=ratio*int(len(array)):
			A=0
		if T<=ratio*int(len(array)):
			T=0
		if C<=ratio*int(len(array)):
			C=0
		if G<=ratio*int(len(array)):
			G=0
		if D<=ratio*int(len(array)):
			D=0

		if A>ratio*int(len(array)):
			NuType+=1
		if T>ratio*int(len(array)):
			NuType+=1
		if C>ratio*int(len(array)):
			NuType+=1
		if G>ratio*int(len(array)):
			NuType+=1
		if D>ratio*int(len(array)):
			NuType+=1
		
		if NuType==1:
			if A==max(A,T,C,G,D):
				consensus.append("A")
			if T==max(A,T,C,G,D):
				consensus.append("T")
			if C==max(A,T,C,G,D):
				consensus.append("C")
			if G==max(A,T,C,G,D):
				consensus.append("G")
			if D==max(A,T,C,G,D):
				consensus.append("-")

		if NuType==2:
			if D==0:
				if A!=0 and G!=0 and C==0 and T==0:
					consensus.append('R')
				if A!=0 and G==0 and C!=0 and T==0:
					consensus.append('M')
				if A!=0 and G==0 and C==0 and T!=0:
					consensus.append('W')
				if A==0 and G!=0 and C!=0 and T==0:
					consensus.append('S')
				if A==0 and G!=0 and C==0 and T!=0:
					consensus.append('K')
				if A==0 and G==0 and C!=0 and T!=0:
					consensus.append('Y')
				
			if D!=0:
				if A!=0 and G==0 and C==0 and T==0:
					consensus.append('a')
				if A==0 and G!=0 and C==0 and T==0:
					consensus.append('g')
				if A==0 and G==0 and C!=0 and T==0:
					consensus.append('c')
				if A==0 and G==0 and C==0 and T!=0:
					consensus.append('t')
	return consensus

def GCcalc(record):
	countMax=0;countMin=0
	for i in range(len(record)):
		if record[i] in {'C','G','S','K','Y','R'}:
			countMax+=1
		if record[i] in {'C','G'}:
			countMin+=1
		GCMean=float((countMin+countMax)/(2*len(record)))
	return ('{:.1%}'.format(GCMean))

def Tmcalc(array):
	pmlist=[]
	for n in range(len(array)):
		if set(array[n])=={'-'}:
			my_tm='NA'
		else:
			my_seq=Seq(''.join(array[n])).upper()
			p1=200;p2=50;p3=2;p4=1;p5=1;preset=mt.DNA_NN4
			my_tm=mt.Tm_NN(my_seq,nn_table=preset,dnac1=p1,Na=p2,Mg=p3,dNTPs=p4,saltcorr=p5)
		pmlist.append(my_tm)
	return (pmlist)


def FindSeq(spec_array,align_array,minsize,maxsize,SpecArrayLength,FullLength,FullArrayLength,degNo,Cutoff):
	seq=[];seqlist=[];seqloc=[];seqGC=[];seqTm=[]
	for windowsize in range(minsize,maxsize+1):
		#Count the nucleotide frequences in the species array
		for n in range(FullLength-windowsize+1):
			CountType=0;degenerate=0
			for nt in range(n,n+windowsize):
				A=0;T=0;C=0;G=0;D=0;NuType=0
				for row in range(SpecArrayLength):
					if spec_array[row][nt]in{'a','A'}:
						A+=1
					elif spec_array[row][nt]in{'t','T'}:
						T+=1
					elif spec_array[row][nt]in{'c','C'}:
						C+=1
					elif spec_array[row][nt]in{'g','G'}:
						G+=1
					else:
						D+=1

				if A<=Cutoff*SpecArrayLength:
					A=0
				if T<=Cutoff*SpecArrayLength:
					T=0
				if C<=Cutoff*SpecArrayLength:
					C=0
				if G<=Cutoff*SpecArrayLength:
					G=0
				if D<=Cutoff*SpecArrayLength:
					D=0

				#Count the types of nucleotides appear in the species array
				if A>Cutoff*SpecArrayLength:
					NuType+=1
				if T>Cutoff*SpecArrayLength:
					NuType+=1
				if C>Cutoff*SpecArrayLength:
					NuType+=1
				if G>Cutoff*SpecArrayLength:
					NuType+=1
				if D>Cutoff*SpecArrayLength:
					NuType+=1
				#Count the conserved location in the species array
				if D==0:
					if NuType<3:
						CountType+=1
					if NuType==2:
						degenerate+=1
				elif D!=0:
					if NuType==1:
						CountType+=1
					else:
						CountType=0
		
			#Use function to generate the consensus sequence of the sliding window in species array
			consensus=[]	
			if CountType==windowsize and degenerate<=degNo:
				consensus=ConsensusSeq(spec_array[:,n:n+windowsize],Cutoff)
				#print(consensus)
				klist=[]
				#To compare this species consensus sequence with all the other aligned sequences 
				#in order to decide if this sequence is specific enough or highly conserved
				for row in range(FullArrayLength):
					k=0;count=0;fullmatch=0
					#Degenerate nucleotides are included
					#Compare nucleotides of the consensus with other sequences and calculate a total number 
					for nt in range(n,n+windowsize):
						if align_array[row][nt]in{'a','A'} and consensus[nt-n] in {'A','R','W','M'}:
							k+=1
						if align_array[row][nt]in{'t','T'} and consensus[nt-n] in {'T','Y','W','K'}:
							k+=1
						if align_array[row][nt]in{'c','C'} and consensus[nt-n] in {'C','Y','S','M'}:
							k+=1
						if align_array[row][nt]in{'g','G'} and consensus[nt-n] in {'G','R','S','K'}:
							k+=1
						if align_array[row][nt]=='-' and consensus[nt-n]=='-':
							k+=1
					klist.append(k)
				#print(klist)
				#Compare the calculation with the criteria set by user
				#Once the calculation matches the criteria, the sequence will be stored
				for item in klist:
					if item==windowsize:
						fullmatch+=1

				if count<=SpecArrayLength or fullmatch==FullArrayLength:
					if set(Tmcalc(spec_array[:,n:n+windowsize]))=={'NA'}:
						seqTm.append('NA')
					else:	
						seqTm.append("%.1f"%(mean(item for item in (Tmcalc(spec_array[:,n:n+windowsize])) if item!='NA')))
					
					seqlist.append(Seq(''.join(consensus)))
					seqGC.append(GCcalc(consensus))
					seqloc.append(n+1)			

	return seqlist,seqloc,seqGC,seqTm
