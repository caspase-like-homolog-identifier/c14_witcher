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


# In[7]:


dyads = prosite.get_dyad()


# In[8]:


dyads


# In[9]:


deathdomains = FindDeathDomains(inseq, "/opt/DB_REF/Pfam/Ig*hmm")


# In[10]:


deathdomains.deathdomains_iter()


# In[11]:


#AVAILABLE FEATURES FOR DEATH DOMAINS
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
dd = deathdomains.DeathDomains('evalue')


# In[17]:


dd


# In[18]:


c14 = C14Subunits(stockholm, dyads)


# In[19]:


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




