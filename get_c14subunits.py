#!/usr/bin/env python
from Bio import pairwise2
from Bio import AlignIO
import operator
import collections
import warnings



class C14Suunitst(object):

    def __init__(self, msa, cys_dict):
        
        self.msa = msa
        self.cys_dict = cys_dict

        
    def cys_location(self):

        stop_lists = []
        for seq in self.msa:
            cys_active_site =  self.cys_dict.get(seq.id, False)
            if cys_active_site:
                #print(cys_active_site)
                alignments = pairwise2.align.localms(seq.seq, cys_active_site, 5, -4, -2, -1, one_alignment_only = True)
                #print(pairwise2.format_alignment(*alignments[0]))
                stop_lists.append(alignments[0][4])
        stop_counter = collections.Counter(stop_lists)
        cys_end = max(dict(stop_counter).items(), key=operator.itemgetter(1))
        counts_totals = sum(stop_counter.values())
        print(stop_counter)
        print("Cysteine position:{:>3}\nConfidence: {:>10.0%}\n".format(cys_end[0],cys_end[1]/float(counts_totals)))

        return cys_end[0]

    def msa_slice(self, start, stop):

        return self.msa[:,start:stop]
        

        
if __name__ == '__main__':
    
    msa =  AlignIO.read("MCA.fasta","fasta")
    c14_data = open("c14.data").read().splitlines()    
    data = dict([ (line.split()[0],line.split()[2]) for line in c14_data if len(line.split()) > 2 ])
    c14 = C14Suunitst(msa, data)
    stop = c14.cys_location()
    print(c14.msa_slice(stop, stop+10))
    
