#!/usr/bin/env python
from Bio.ExPASy import Prosite
import argparse
import pattern
import re



class PrositePattern(object):

    def __init__(self):

        """instatiates prosite pattern search"""
        
        #https://prosite.expasy.org/PS01121
        self.HIS_pattern = "H(0,1)-x(2,4)-[SCQ](0,1)-x(2)-{A}-x-[LIVMF](0,3)-[HFYST](0,2)-[GST](1,2)-H-G"
        #https://prosite.expasy.org/PS01122
        self.CYS_pattern = "[STAVK]-[HTRKP]-[KIMLV]-[LIVMFHT]-[LIVMFYSCA]-[SLIVMF](2)-[DQP]-[TCSAF]-C-[RHQG]-[SGE]"

        self.regex_HIS = re.compile(pattern.prosite_to_re(self.HIS_pattern))
        self.regex_CYS = re.compile(pattern.prosite_to_re(self.CYS_pattern))

        
    def Caspase_HIS(self, seq):
        """find p20 HIS residue"""
        return self.regex_HIS.findall(seq)

    
    def Caspase_CYS(self, seq):
        """find p20 HIS residue"""
        return self.regex_CYS.findall(seq)

    
if  __name__ ==  '__main__':
    from Bio import SeqIO
    parser = argparse.ArgumentParser("Find P20 Cystine and Histidine dyads")
    parser.add_argument('sequences',help ="sequences file", type=argparse.FileType('r'))
    parser.add_argument('-f','--format', default = "fasta" ,help ="sequences format")
    args = parser.parse_args()

    prosite = PrositePattern()
    records = SeqIO.parse(args.sequences, args.format)
    
    j = 0
    for i,rec in enumerate(records, 1):
        cys = "".join(prosite.Caspase_CYS(str(rec.seq.ungap('-'))))
        his = "".join(prosite.Caspase_HIS(str(rec.seq.ungap('-'))))
        print("{:<20}{:<35}{:<50}".format(rec.id,his,cys))
    
