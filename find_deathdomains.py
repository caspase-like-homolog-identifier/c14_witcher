#!/usr/bin/env python

from run_hmmer import RunHmmer
import argparse
import glob


class FindDeathDomains(RunHmmer):

     def __init__(self, seqfile, dd_hmm_path, *hmmersearch_args):
         
         self.dd_hmm_paths = glob.glob(dd_hmm_path)
         super().__init__("hmmsearch", None, seqfile, None, None, *args)
         
         



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('seqfile', action='store', type=str)
    parser.add_argument('-g','--hmm_glob', default="/opt/DB_REF/Pfam/Ig*hmm")
    args = parser.parse_args()    
    dd = FindDeathDomains(args.seqfile, args.hmm_glob)
    
