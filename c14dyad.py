#!/usr/bin/env python
from Bio.ExPASy import Prosite
from Bio import SeqIO
import pattern
import re

# {PS01121; CASPASE_HIS}
# {PS01122; CASPASE_CYS}
# {PS50207; CASPASE_P10}
# {PS50208; CASPASE_P20}

class PrositePattern(object):

    def __init__(self, prosite_db = "prosite.dat"):

        #H-x(2,4)-[SC]-x(2)-{A}-x-[LIVMF](2)-[ST]-H-Ge
        #K-P-K-[LIVMF]-[LIVMFY]-[LIVMF](2)-[QP]-[AF]-C-[RQG]-[GE].
        self.CYS.pattern = "K-P-K-[LIVMF]-[LIVMFY]-[LIVMF](2)-[QP]-[AF]-C-[RQG]-[GE]."
        self.CYS.pattern = "[STAVK]-[HTRKP]-[KIMLV]-[LIVMFHT]-[LIVMFYSCA]-[SLIVMF](2)-[DQP]-[TCSAF]-C-[RHQG]-[SGE]"
    
        self.regex_HIS = re.compile(pattern.prosite_to_re(self.HIS.pattern))
        self.regex_CYS = re.compile(pattern.prosite_to_re(self.CYS.pattern))
                    
    def Caspase_HIS(self, seq):
        return self.regex_HIS.findall(seq)


    def Caspase_CYS(self, seq):
        return self.regex_CYS.findall(seq)



if  __name__ ==  '__main__':
    prosite = PrositePattern()
    records = SeqIO.parse("MCA.fasta","fasta")
    #records = SeqIO.parse("PF00656_seed_align.txt","fasta")
    j = 0
    for i,rec in enumerate(records, 1):
        find = prosite.Caspase_CYS(str(rec.seq.ungap('-')))
        
        print(i, j, rec.id, find)
        
        if find:
            j += 1
            
    print("\n\n{:.2%} ({}/{})\n".format(j/float(i),j,i))
            

        
    #alignments = next(records)
    #for aln in alignments:
        #print(aln)
        #seq = str(aln.seq.ungap("-").upper())
        #print(prosite.Caspase_CYS(seq))
    #prosite.
    #https://docs.python.org/3/library/re.html
    #print(pattern.prosite_to_re(self.HIS))
    #patt = "H-x(2,4)-[SC]-x(2)-{A}-x-[LIVMF](2)-[ST]-H-G"
    #patt = "H.{2,4}[SC].{2}[^A].[LIVMF]{2}[ST]HG"
      
