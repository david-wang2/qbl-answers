#!/bin/bash

cat week4_quer.fa seqdump.txt > seqdump2.txt
transeq seqdump2.txt
mafft query.pep > query.dat
