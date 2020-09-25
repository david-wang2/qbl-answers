#!/bin/bash

for NUM in 09 11 23 24 27 31 35 39 62 63
do
bwa mem -R "@RG\tID:A01_${NUM}\tSM:A01_${NUM}" sacCer3.fa A01_${NUM}.fastq > aln-A01_${NUM}.sam
done
