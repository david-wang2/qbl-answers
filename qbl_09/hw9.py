#!/usr/bin/env python2
from __future__ import division
import hifive
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig

# Get data
hic=hifive.HiC("filtered_1.dat",'r')
chr13=hic.cis_heatmap('chr13',1000000,datatype='fend',arraytype='full',diagonalincluded=True)
enrichment=(chr13[:,:,0]+1)/(chr13[:,:,1]+1)
log_enrichment=np.log(enrichment)

# Create heatmap of the log of enrichment scores
fig,ax=plt.subplots(figsize=(14,10))
ax.set_title("Heatmap of Enrichment Scores for Chr13",fontsize=20)
ax=sns.heatmap(log_enrichment)
plt.savefig("chr13_heatmap.png")

# Compartment Analysis 
Comp = hifive.hic_domains.Compartment(hic, 100000, chroms=['chr13'], out_fname='tmp.hdf5')
Comp.write_eigen_scores('hic_comp.bed')
X = Comp.positions['chr13']
Y = Comp.eigenv['chr13']

fig,ax=plt.subplots(figsize=(14,10))
plt.rcParams['xtick.labelsize']=20
plt.rcParams['ytick.labelsize']=20
plt.scatter(X[:,0],Y)
ax.set_title("Comparment Scores",fontsize=20)
ax.set_xlabel("Position (Mb)",fontsize=20)
ax.set_ylabel("Eigenscores",fontsize=20)
plt.savefig("Comparment_scores.png")

# separate genes and create a violin plot
genes = {}
with open('data/WT_fpkm.bed','r') as f:
    for line in f:
        data = line.split('\t')
        chromosome = data[0]
        name=data[3]
        start = int(data[1])
        end = int(data[2])
        expression = float(data[4].split('\n')[0])
        genes[name]=(name,start,end,expression)

# Compute total length of a gene as there are multiple copies
expressed = {}
with open('comp_positive_genes.bed','r') as f:
    for line in f:
        data =line.split('\t')
	gene=data[3]
	start=int(data[1])
	end=int(data[2])
        print(gene)
	if gene in expressed.keys():
            expressed[gene]+=end-start
	else:
            expressed[gene]=end-start

expressed_genes = []
for gene in expressed:
    if gene in genes:
        size=expressed[gene]/(genes[gene][2]-genes[gene][1])
    # crude filter for size
        if size>0.5: 
            expressed_genes.append(genes[gene][3])

repressed = {}
with open('comp_negative_genes.bed','r') as f:
    for line in f:
        data =line.split('\t')
        gene=data[3]
        start=int(data[1])
        end=int(data[2])
        if gene in repressed.keys():
            repressed[gene]+=end-start
        else:
            repressed[gene]=end-start

repressed_genes = []
for gene in repressed:
    if gene in genes:
        size=repressed[gene]/(genes[gene][2]-genes[gene][1])
        # crude filter for size
        if size>0.5:
            repressed_genes.append(genes[gene][3])

fig, ax = plt.subplots()
ax.violinplot(dataset=[expressed_genes, repressed_genes],showmedians=True)
ax.set_xticks([1,2])
ax.set_xticklabels(["Compartment A","Compartment B"],fontsize=20)
ax.set_ylabel("Expression",fontsize=20)
ax.set_title("Gene Expression in Compartments A and B",fontsize=20)
plt.savefig("violin.png")
