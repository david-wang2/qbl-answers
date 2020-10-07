import sys
from fasta_parser import FASTAReader
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# read nc_align.dat, seqdump.txt

# Question 3
nucSeq = FASTAReader(open("seqdump2.txt","r"))
protAlign = FASTAReader(open("query.dat",'r'))

outStr = ""
for i in range(len(nucSeq)):
    nucAlign = nucSeq[i][0] + '\n'
    nucCount = 0
    for aa in protAlign[i][1]:
        if (aa == '-'):
            nucAlign += '---'
        else:
            nucAlign += nucSeq[i][1][nucCount:nucCount+3]
            nucCount += 3
    nucAlign += '\n'
    outStr += nucAlign

with open("nuc_align.dat",'w') as f:
    f.write(outStr)

# Question 4
codon = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

# Read the output from question 3 and don't include header
count = 1
lstNuc = []
with open("nuc_align.dat",'r') as f:
    for line in f:
        # if test to test for header
        if (count % 2 == 0):
            lstNuc.append(line.split('\n')[0])
        count += 1

# Parse the alignment file and find synonymous and nonsynonymouus mutations
blast_seq1 = lstNuc[0]
seqLen = int(len(blast_seq1)/3)
dN = 0
lst_dN = []
dS = 0
lst_dS = []

count = -1
for i in range(seqLen):
    if (blast_seq1[i*3:(i+1)*3] == '---'):
        pass
    else:
        count += 1
        s = 0
        n = 0
        for j in range(1,len(lstNuc)):
            nuc1 = blast_seq1[i*3:(i+1)*3]
            nuc2 = lstNuc[j][i*3:(i+1)*3]
            if (nuc1 == '---') or (nuc2 == '---') or ('R' in nuc1) or ('R' in nuc2) or ('Y' in nuc1) or ('Y' in nuc2) or ('S' in nuc1) or ('S' in nuc2) or ('N' in nuc1) or ('N' in nuc2) or ('W' in nuc1) or ('W' in nuc2):
                break
            elif (nuc1 != nuc2) and (codon[nuc1] == codon[nuc2]):
                s += 1
                dS += 1
            elif (nuc1 != nuc2) and (codon[nuc1] != codon[nuc2]):
                n += 1
                dN += 1
        lst_dS.append(s)
        lst_dN.append(n)

# calculate z-scores for each codon
lst_D = np.subtract(np.array(lst_dN),np.array(lst_dS))
ste = np.std(lst_D)/np.sqrt(len(lst_D))
zScores = []
for i in range(len(lst_D)):
    zScores.append(lst_D[i]/ste)
pValue = stats.norm.sf(abs(np.array(zScores)))*2

# create ratio for plotting and mark as significant
poi = pValue < 0.05
sigP_x = []
sigP_y = []
lstCodon = []
outDat = []
for i in range(len(lst_dS)):
    if (lst_dS[i] != 0) and (lst_dN[i] != 0):
        lstCodon.append(i+1)
        outDat.append(np.log2(lst_dN[i]/lst_dS[i]))
        if poi[i]:
            sigP_x.append(i+1)
            sigP_y.append(np.log2(lst_dN[i]/lst_dS[i]))

# Create scatterplot of log ratio vs codon with significant and nonsignificant points
plt.figure(figsize=(14,10))
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.scatter(lstCodon,outDat, color = "purple", label = "Nonsignificant")
plt.scatter(sigP_x,sigP_y,color = "green", label = "Significant")
plt.legend(loc='upper right',fontsize=20)
plt.xlabel('codon', fontsize = 20)
plt.ylabel(r'log$_2$($\dfrac{dN}{dS}$)', fontsize = 20)
plt.title('Log ratio of Synonymous to Nonsynonymous mutations', fontsize = 24)
plt.savefig('mutation_ratio.png')
plt.close('all')
