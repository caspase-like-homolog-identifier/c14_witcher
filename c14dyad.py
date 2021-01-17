#!/usr/bin/env python
from Bio import SeqIO
import pandas as pd
import argparse
import pattern
import re



class PrositePattern(object):

    def __init__(self, alignment_fname, msa_format="stockholm", gap_char = "-"):

        """instatiates prosite pattern search"""
        
        #https://prosite.expasy.org/PS01121
        self.HIS_pattern = "H(0,1)-x(2,4)-[SCQ](0,1)-x(2)-{A}-x-[LIVMF](0,3)-[HFYST](0,2)-[GST](1,2)-H-G"
        #https://prosite.expasy.org/PS01122
        self.CYS_pattern = "[STAVK]-[HTRKP]-[KIMLV]-[LIVMFHT]-[LIVMFYSCA]-[SLIVMF](2)-[DQP]-[TCSAF]-C-[RHQG]-[SGE]"

        self.regex_HIS = re.compile(pattern.prosite_to_re(self.HIS_pattern))
        self.regex_CYS = re.compile(pattern.prosite_to_re(self.CYS_pattern))
        self.alignment = list(SeqIO.parse(alignment_fname, msa_format))
        self.gap_char = gap_char
        
    def Caspase_HIS(self, seq):
        """find p20 HIS residue"""
        
        return self.regex_HIS.findall(seq)

    
    def Caspase_CYS(self, seq):
        """find p20 HIS residue"""
        
        return self.regex_CYS.findall(seq)

    
    def get_dyad(self):
        """Get a data of the dyad activate sites"""


        columns = ['Seq_ID', 'Caspase_CYS', 'CASPASE_HIS']
        dyads =  pd.DataFrame(columns = columns)
        for i,seq in enumerate(self.alignment, 1):
             cys = "".join(self.Caspase_CYS(str(seq.seq.ungap(self.gap_char).upper())))
             his = "".join(self.Caspase_HIS(str(seq.seq.ungap(self.gap_char).upper())))
             dyads = dyads.append(pd.Series(index = columns, name = seq.id, data=[seq.id, cys, his]))
             
        return dyads
            
            
    
if  __name__ ==  '__main__':
    from Bio import SeqIO
    parser = argparse.ArgumentParser("Find P20 Cystine and Histidine dyads")
    parser.add_argument('alignment',help ="sequences file", type=argparse.FileType('r'))
    parser.add_argument('-f','--format', default = "stockholm" ,help ="sequences format")
    parser.add_argument('-g','--gap_char', default = "-" ,help ="gap character in alignment")
    args = parser.parse_args()
    prosite = PrositePattern(args.alignment, args.format, args.gap_char)
    prosite.get_dyad()
        
