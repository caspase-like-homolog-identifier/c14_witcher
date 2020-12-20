#!/usr/bin/env python
from Bio.ExPASy import Prosite
from Bio import AlignIO
import pattern
import re


# {PS01121; CASPASE_HIS}
# {PS01122; CASPASE_CYS}
# {PS50207; CASPASE_P10}
# {PS50208; CASPASE_P20}

class PrositePattern(object):

    def __init__(self, prosite_db = "prosite.dat"):

        prosite_handle = open(prosite_db)
        prosite_records = Prosite.parse(prosite_handle)
        prosite_dict = dict([ (rec.accession, rec) for rec in prosite_records ])

        self.HIS = prosite_dict["PS01121"]
        self.CYS = prosite_dict["PS01122"]
        
        self.regex_HIS = re.compile(pattern.prosite_to_re(self.HIS.pattern))
        self.regex_CYS = re.compile(pattern.prosite_to_re(self.CYS.pattern))
                

        
    def Caspase_HIS(self, seq):
        return self.regex_HIS.findall(seq)


    
    def Caspase_CYS(self, seq):
        return self.regex_CYS.findall(seq)



if  __name__ ==  '__main__':
    prosite = PrositePattern()
    records = AlignIO.parse("in_file_trim.ali","stockholm")
    alignments = next(records)
    for aln in alignments:
        seq = str(aln.seq.ungap("-").upper())
        print(prosite.Caspase_CYS(seq))
    #prosite.
    #https://docs.python.org/3/library/re.html
    #print(pattern.prosite_to_re(self.HIS))
    #patt = "H-x(2,4)-[SC]-x(2)-{A}-x-[LIVMF](2)-[ST]-H-G"
    #patt = "H.{2,4}[SC].{2}[^A].[LIVMF]{2}[ST]HG"
      
