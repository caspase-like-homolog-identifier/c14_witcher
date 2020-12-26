#!/usr/bin/env python
from Bio.Align import AlignInfo
from Bio import pairwise2
from Bio import AlignIO
import pprint
import operator
import collections
import warnings



class C14Subunits(object):

    def __init__(self, msa, dyad_dict, consensus_threshold = 0.3):
        
        self.msa = msa
        self.dyad_dict = dyad_dict 
        align_info = AlignInfo.SummaryInfo(self.msa)        
        self.align_consensus = align_info.gap_consensus(threshold = consensus_threshold)

        cys_lists = []
        his_lists = []
        for seq in self.msa:
            his_site, cys_site = self.dyad_dict.get(seq.id, ("-","-",))
            if cys_site != '-':
                  cys_align = pairwise2.align.localms(seq.seq, cys_site, 5, -4, -2, -1, one_alignment_only = True)
                  cys_lists.append(cys_align[0][4])
                  #print(seq.id)
                  #print(pairwise2.format_alignment(*cys_align[0]))
            if his_site != '-':
                  his_align = pairwise2.align.localms(seq.seq, his_site, 5, -4, -2, -1, one_alignment_only = True)
                  his_lists.append(his_align[0][4])
        
        self.his_end = self.get_conf("Histidine", his_lists)
        self.cys_end = self.get_conf("Cysteine", cys_lists)
        self.col_end = False
        
    def get_conf(self, site, pos_lists):
                  
          if not pos_lists:
              return None
          stop_counter = collections.Counter(pos_lists)
          end = max(dict(stop_counter).items(), key=operator.itemgetter(1))
          counts_totals = sum(stop_counter.values())
          print(stop_counter)
          print("{} position: {}\nConfidence: {:.0%}\n".format(site, end[0],end[1]/float(counts_totals)))
    
          return end[0]

            
    def get_p20(self):

        if not (self.cys_end):
            print(self.cys_end)

        #Find the find aspartic acid (D) after the catalytic active site (C)
        self.col_end = min([ pos for pos,aa in enumerate(str(self.align_consensus),1) if (aa == 'D' and pos > self.cys_end)])

        #print("**", self.col_end)
        return self.msa[:,0:self.col_end]

    
    def get_p10(self):
        
        if not self.col_end:
            self.get_p20()
        
        score_tr = {'.': -1,'*':10 }
        column_pp = { i: int(score_tr.get(score,score)) for i,score in enumerate(self.msa.column_annotations['posterior_probability'],1) }
        #pprint.pprint(column_pp)
        
       
if __name__ == '__main__':

    data_fname = "MCA2.data"
    msa_fname = "MCA2.sto"
    dyad_dict = {}    
    with open(data_fname) as fp:
        for line in fp:
             line = line.split()
             data = line[1:]
             if data:
                dyad_dict[line[0]] = data

    msa =  AlignIO.read(msa_fname, "stockholm")
    c14 = C14Subunits(msa, dyad_dict)
    print(c14.get_p10())
    #stop = c14.cys_location()
    #print(c14.msa_slice(stop, stop+10))
    
