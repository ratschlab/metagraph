#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


path = "/cluster/home/mmarzett/superkmer_stats_color.txt"


# In[3]:


file = open(path, "r")
colors_array = np.loadtxt(path, dtype=int, skiprows=1)


# In[4]:


colors_array


# In[5]:


plt.hist(colors_array,bins=np.arange(min(colors_array), max(colors_array)+1, step=1))
plt.title("#color changes per superkmer")
plt.xlabel("#color changes")
plt.savefig("colors_hist_10000")
plt.show()


# In[ ]:




