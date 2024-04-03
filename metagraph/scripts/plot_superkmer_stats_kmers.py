#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


path = "/cluster/home/mmarzett/superkmer_stats_kmers.txt"


# In[10]:


file = open(path, "r")
kmers_array = np.loadtxt(path, dtype=int, skiprows=1)


# In[11]:


kmers_array


# In[12]:


plt.hist(kmers_array,bins=np.arange(min(kmers_array), max(kmers_array)+1, step=1))
plt.title("#kmers per superkmer")
plt.xlabel("#kmers")
plt.savefig("kmers_hist")
plt.show()


# In[ ]:




