# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from find_deathdomains import FindDeathDomains
from get_c14subunits import C14Subunits
from decisiontree import DecisionTree
from c14dyad import PrositePattern
from run_hmmer import RunHmmer
import pandas as pd
import collections

inseq = "MCAx.fasta"
stockholm = "MCAx.sto"

hmmsearch = RunHmmer('hmmsearch',
                     'Peptidase_C14.hmm', 
                     inseq,
                     stockholm,
                     "MCA.tbloutx",
                     "--cpu 4", 
                     "--incE 0.00001")

out, err = hmmsearch()

prosite = PrositePattern(stockholm)

dyads = prosite.get_dyad()

deathdomains = FindDeathDomains(inseq, "/opt/DB_REF/Pfam/Ig*hmm")



deathdomains.deathdomains_iter()

deathdomains.print_deathdomains()

dd_dict = deathdomains.DeathDomains

dd_dict

dyads

FindDeathDomains.dumb_merge(dyads, dd_dict)

dd

c14 = C14Subunits(stockholm, dyads)

p10 = c14.get_p10()

print(p10)

p20 = c14.get_p20()

print(p20)

linker = c14.get_linker()

print(linker)

stats = c14.get_stats()

stats

dt = DecisionTree(stats)

dt.classify()


