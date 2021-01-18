#!/usr/bin/env python

from run_hmmer import RunHmmer
from Bio import SearchIO
import collections
import random
import tempfile
import argparse
import pprint
import glob
import sys

class FindDeathDomains(RunHmmer):

     def __init__(self, seqfile, dd_hmm_path, *hmmersearch_args):

         """ Subclass the Hmmer commandline wrapper """ 
          
         self.dd_hmm_paths = glob.glob(dd_hmm_path)
         super().__init__("hmmsearch", None, seqfile, None, None, *hmmersearch_args)
         self.deathdomain_hits = {}
         self.dd_dict = None

     def deathdomains_iter(self):

          """ iterate over the deathdomains """
          dd_names  = set()
          for hmm_file in self.dd_hmm_paths:
               self.hmmfile = hmm_file
               tmp1, tmp2 = [ tempfile.NamedTemporaryFile(delete=False) for _ in range(2) ]
               self.align_out = tmp1.name
               self.domtblout = tmp2.name
               std, stderr = self()
               deathdomain = self.has_deathdomain(self.domtblout)

               if deathdomain:
                    self.deathdomain_hits[deathdomain[0].id] = deathdomain[0].hits 
                    dd_names.add(deathdomain[0].id)
               
                   
     def has_deathdomain(self, domtab):

          return list(SearchIO.parse(domtab, "hmmsearch3-domtab"))



     def print_deathdomains(self):
     
         if not self.dd_dict:
               return None
          
         for k,v in self.dd_dict.items():
             print("{}\t{}".format(k,"\t".join(v)))
          
          
     
     def DeathDomains(self, feature):
         """Property to view the death domains.Start analysis if not done already"""
         # _id
         # _id_alt
         # _query_id
         # _description
         # _description_alt
         # _query_description
         # attributes
         # dbxrefs
         # _items
         # accession
         # seq_len
         # evalue
         # bitscore
         # bias
         
         if not self.deathdomain_hits:
              self.deathdomains_iter()
         #create dict using seq.ids as keys and empty lists as values
         self.dd_dict = collections.defaultdict(list)
         for dd in self.deathdomain_hits:
             #print(dd)
             for hit in self.deathdomain_hits[dd]:
                 self.dd_dict[hit.id].append(vars(hit)[feature])

         pprint.pprint(pd.DataFrame(self.dd_dict))
         
         #return self.dd_dict                    

     
     def c14merge(c14seq_dict1, c14seq_dict2):
          
         merged = c14seq_dict1.copy()
         
         for key in merged:
             try:
                merged[key] = merged[key] + c14seq_dict2.get(key,[])
             except KeyError:
                merged[key] = c14seq_dict2.get(key,[])

         #Make sure we did not miss any keys in dict2
         merged.update(c14seq_dict2)
                  
         return merged
               
 
     
     
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('seqfile', action='store', type=str)
    parser.add_argument('-g','--hmm_glob', default="/opt/DB_REF/Pfam/Ig*hmm")
    args = parser.parse_args()    
    dd = FindDeathDomains(args.seqfile, args.hmm_glob)
    dd.deathdomains_iter()
    print("\n\n\n\n")
    
    print(dd.DeathDomains('evalue'))
