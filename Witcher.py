#!/usr/bin/env python
# coding: utf-8

# In[1]:


from find_deathdomains import FindDeathDomains
from get_c14subunits import C14Subunits
from decisiontree import DecisionTree
from c14dyad import PrositePattern
from run_hmmer import RunHmmer
import pandas as pd
import collections


# In[2]:


df = pd.DataFrame(columns=['a','b','c'])


# In[3]:


inseq = "MCAx.fasta"
stockholm = "MCAx.sto"


# In[4]:


hmmsearch = RunHmmer('hmmsearch',
                     'Peptidase_C14.hmm', 
                     inseq,
                     stockholm,
                     "MCA.tbloutx",
                     "--cpu 4", 
                     "--incE 0.00001")


# In[5]:


out, err = hmmsearch()


# In[6]:


prosite = PrositePattern(stockholm)


# In[15]:


dyads = prosite.get_dyad()


# In[16]:


dyads


# In[17]:


deathdomains = FindDeathDomains(inseq, "/opt/DB_REF/Pfam/Ig*hmm")


# In[18]:


deathdomains.deathdomains_iter()


# In[26]:


deathdomains.print_deathdomains()


# In[27]:


dd_dict = deathdomains.DeathDomains


# In[29]:


vars(str)


# In[28]:


dd_dict


# In[14]:


FindDeathDomains.dumb_merge(dyads, dd_dict)


# In[ ]:


dd


# In[ ]:


c14 = C14Subunits(stockholm, dyads)


# In[ ]:


p10 = c14.get_p10()


# In[ ]:


print(p10)


# In[ ]:


p20 = c14.get_p20()


# In[ ]:


print(p20)


# In[ ]:


linker = c14.get_linker()


# In[ ]:


print(linker)


# In[ ]:


stats = c14.get_stats()


# In[ ]:


stats


# In[ ]:


dt = DecisionTree(stats)


# In[ ]:


dt.classify()


# In[ ]:




