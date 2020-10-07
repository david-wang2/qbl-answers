#!usr/bin/env python2

import sys

def FASTAReader(file):
    line = file.readline()
    
    assert line.startswith('>'), 'Not a fasta file'
    
    seq_id = line[1:].strip('\n\r')
    
    sequence = []
    sequences = []
    line = file.readline()
    while line: # while haven't reached the end of the file
        if line.startswith('>'):
            sequences.append((seq_id,''.join(sequence)))
            seq_id = line[1:].rstrip('\n\r')
            sequence = []
        else:
            sequence.append(line.strip())
        line = file.readline()
        
    sequences.append((seq_id,''.join(sequence)))
    seq_id = line[1:].rstrip('\n\r')
    return sequences

if __name__ == "__main__":
    seqs = FASTAReader(sys.stdin)
    for seq_id,sequence in seqs:
        print("{}: {}".format(seq_id,sequence),file=sys.stdout)
