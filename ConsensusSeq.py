from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

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
