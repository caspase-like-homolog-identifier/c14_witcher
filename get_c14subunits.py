#!/usr/bin/env python
from Bio import pairwise2
from Bio import AlignIO
import collections
import warnings




class C14Suunitst(object):

    def __init__(self, align_sequence, sub_seq):
        
        self.align_sequence = align_sequence
        self.sub_seq = sub_seq
        
    def local_align(self):

        alignments = pairwise2.align.localms(self.align_sequence, self.sub_seq, 5, -4, -2, -1, one_alignment_only = True)
        #print(pairwise2.format_alignment(*alignments[0]))
        
        return alignments[0][4]
        
        
if __name__ == '__main__':
    msa =  AlignIO.read("MCA.fasta","fasta")
    c14_data = open("c14.data").read().splitlines()
    data = dict([ (line.split()[0],line.split()[1]) for line in c14_data if len(line.split()) > 1 ])
    i = 0
    stop_lists = []
    for seq in msa:
        dyad = data.get(seq.id,False)
        if dyad:
            #print(seq.seq)
            c14 = C14Suunitst(str(seq.seq).replace("-",'.'), dyad)
            stop = c14.local_align()
            stop_lists.append(stop)
        i += 1
        if i == 5:
            break
        
    print(collections.Counter(stop_lists))

        
