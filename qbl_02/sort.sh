#!/bin/bash

for NUM in 09 11 23 24 27 31 35 39 62 63
do
samtools sort aln-A01_${NUM}.sam -o A01_${NUM}.bam -O bam
done
