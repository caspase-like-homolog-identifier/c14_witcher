#!/usr/bin/env python
from Bio.Align import AlignInfo
from Bio import pairwise2
from Bio import AlignIO
import pandas as pd
import numpy as np
import pprint
import operator
import collections
import warnings



class C14Subunits(object):

    def __init__(self, msa_fname, dyad_df, msa_format="stockholm", consensus_threshold = 0.3):

        """Get alignment and dataframe correspondind dyads and analyse for p10 and p20 subunit as well a linker"""
        
        self.msa = msa = AlignIO.read(msa_fname, msa_format)
        dyad_df.index = list(map(lambda index: index.split('/')[0], dyad_df.Seq_ID))
        align_info = AlignInfo.SummaryInfo(self.msa)        
        self.align_consensus = align_info.dumb_consensus(threshold = consensus_threshold)
        score_tr = {'.': -1,'*':10 }
        self.column_pp = { i: int(score_tr.get(score,score)) for i,score in enumerate(self.msa.column_annotations['posterior_probability'],1) }
        self.msa_len = self.msa.get_alignment_length()
        self.p10 = None 
        self.p20 = None
        cys_lists = []
        his_lists = []
        for seq in self.msa:
            if seq.id in dyad_df.index:         
               his_site, cys_site = dyad_df.loc[seq_id,["Caspase_CYS", "CASPASE_HIS"]]
               if cys_site != np.nan:
                   print(seq.id, str(cys_site))
                   cys_align = pairwise2.align.localms(seq.seq, str(cys_site), 5, -4, -2, -1, one_alignment_only = True)
                   print(cys_align)
                   aa_seq, stop = cys_align[0][0:5:4]
                   self.get_position(aa_seq, stop, "C")
                   cys_lists.append(cys_align[0][4])
               if his_site != np.nan:
                   his_align = pairwise2.align.localms(seq.seq, str(his_site), 5, -4, -2, -1, one_alignment_only = True)
                   his_lists.append(his_align[0][4])

        self.his_end = self.get_conf("Histidine", his_lists)
        self.cys_end = self.get_conf("Cysteine", cys_lists)



    def get_position(self, aa_seq, stop, residue):

        try:
            #print(aa_seq.index(residue, stop))
            print(aa_seq[stop:])
        except ValueError:
            pass
            #print(aa_seq[stop:])
        
        
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
        self.p20 = min([ pos for pos,aa in enumerate(str(self.align_consensus),1) if (aa == 'D' and pos > self.cys_end)])

        #print("**", self.col_end)
        return self.msa[:,0:self.p20]

    
    def get_p10(self, min_score=8):

        """get the start of the p10 subunit """
        
        if not self.p20:
            self.get_p20()
    
        for pos in range(self.p20, self.msa_len):
            pp_score = self.column_pp[pos]
            consensus_aa = self.align_consensus[pos]
            #print(pos, pp_score, consensus_aa)
            if (pp_score >= min_score) and (consensus_aa in "VAI"):
                self.p10 = pos
                return self.msa[:,self.p10:]
            

    def get_linker(self):
        
        """get the interdomain region"""
        
        if not self.p10:
            self.get_p10()
        return self.msa[:,self.p20:self.p10]

       
    def get_length(self, region):
        """ get ungapped sequence length"""
        
        return {seq.id: len(str(seq.seq.ungap('-')).replace('.','')) for seq in region }


    def get_stats(self, fix_ids = False):

        """
            Generate length statistics for all regions 
            Fix ids removes region coordinates from ids and will identify duplicates ids
        """
        
        p10_len = self.get_length(self.get_p10())
        p20_len = self.get_length(self.get_p20())
        linker_len = self.get_length(self.get_linker())
        p10_len_df = pd.DataFrame.from_dict(p10_len,  orient='index', columns = ["p10"])
        p20_len_df = pd.DataFrame.from_dict(p20_len,  orient='index', columns = ["p20"])
        linker_len = pd.DataFrame.from_dict(p20_len,  orient='index', columns = ["linker"])
        
        if not fix_ids:
            return pd.concat([p20_len_df, linker_len, p10_len_df], axis=1)
        
        else:
            length_stats_df = pd.concat([p20_len_df, linker_len, p10_len_df], axis=1)    
            new_index = []
            id_counts =  collections.defaultdict(int)
            for ref in length_stats_df.index:
                ref = ref.split("/")[0]
                id_counts[ref] += 1
                tmp_id =  id_counts[ref]
                if (tmp_id != 1):
                    ref = "_".join(map(str, [ref, tmp_id]))

                new_index.append(ref)
                    
            length_stats_df.index = new_index
            
        return length_stats_df
            
        
                          
if __name__ == '__main__':

    data_fname = "dyad.tsv"
    msa_fname = "MCAx.sto"
    
    dyad_df = pd.read_csv(data_fname, sep ="\t" )
    
        # for line in fp:
        #      line = line.split()
        #      data = line[1:]
        #      if data:
        #         dyad_df[line[0]] = data
                
    c14 = C14Subunits(msa_fname, dyad_df)
    #print(c14.get_stats(fix_ids = True))
    

