#!/bin/bash

# Part 1
bowtie2-build chr19.fa g1e

cd g1e

macs2 callpeak -t CTCF_ER4.bam -c input_ER4.bam --format=BAM -name=ER4
macs2 callpeak -t CTCF_G1E.bam -c input_G1E.bam --format=BAM -name=G1E

bedtools subtract -a ame=ER4_peaks.narrowPeak -b ame=G1E_peaks.narrowPeak > gain.bed
bedtools subtract -a ame=G1E_peaks.narrowPeak -b ame=ER4_peaks.narrowPeak > loss.bed

for feature in promoter intron exon
do
	grep "${feature}" Mus_musculus.GRCm38.94_features.bed.txt>${feature}.bed
done
for ct in G1E ER4
do
	for feature in promoter intron exon
	do
		bedtools intersect -a CTCF_${ct}_peaks.narrowPeak -b ${feature}.bed  >${ct}_${feature}_overlap
	done
done

# Part 2

# 2.1
sort -r -k 9 ER4_peaks.narrowPeak | head -n100 > meme_dis.dat
bedtools getfasta -fi chr19.fa -bed meme_dis.bed > meme_seqs.fa
meme-chip -meme-maxw 20 -oc meme_data meme_seqs.fa

# 2.2
cp motif_databases/JASPAR/JASPAR_CORE_2016.meme ./
tomtom memechip_out JASPAR_CORE_2016.meme
