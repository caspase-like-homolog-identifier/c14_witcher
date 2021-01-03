#!/usr/bin/env python

from run_hmmer import RunHmmer
from Bio import SearchIO
import tempfile
import argparse
import pprint
import glob


class FindDeathDomains(RunHmmer):

     def __init__(self, seqfile, dd_hmm_path, *hmmersearch_args):

         """ Subclass the Hmmer commandline wrapper """ 
          
         self.dd_hmm_paths = glob.glob(dd_hmm_path)
         super().__init__("hmmsearch", None, seqfile, None, None, *hmmersearch_args)
         self.deathdomain_hits = {}
         

     def get_deathdomains(self):

          """ iterate over the deathdomains """

          for hmm_file in self.dd_hmm_paths:
               self.hmmfile = hmm_file
               tmp1, tmp2 = [ tempfile.NamedTemporaryFile(delete=False) for _ in range(2) ]
               self.align_out = tmp1.name
               self.domtblout = tmp2.name
               #print(self.hmmfile)
               std, stderr = self()
               deathdomain = self.has_deathdomain(self.domtblout)
               
               if deathdomain:
                   self.deathdomain_hits[deathdomain[0]] = deathdomain[0].hits 
                   
                    
               

     def has_deathdomain(self, domtab):

          return list(SearchIO.parse(domtab, "hmmsearch3-domtab"))


     def deathdomains(self):

          if not self.deathdomain_hits:
               self.get_deathdomains()
          



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('seqfile', action='store', type=str)
    parser.add_argument('-g','--hmm_glob', default="/opt/DB_REF/Pfam/Ig*hmm")
    args = parser.parse_args()    
    dd = FindDeathDomains(args.seqfile, args.hmm_glob)
    dd.iterater_deathdomains()
