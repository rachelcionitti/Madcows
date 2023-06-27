#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import fits 
import matplotlib.pyplot as plt
import numpy as np 
from scipy.optimize import minimize, rosen, rosen_der
from scipy import interpolate
import scipy
import os 
from scipy.optimize import curve_fit
from scipy.integrate import simps
import math
from astropy.cosmology import WMAP9 as cosmo
from astroquery.mast import Observations


# In[ ]:


#https://mast.stsci.edu/api/v0/_c_a_o_mfields.html
#https://astroquery.readthedocs.io/en/latest/mast/mast.html#id1
#https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
#https://outerspace.stsci.edu/display/MASTDOCS/Search+a+List+of+Targets
#https://www.stsci.edu/hst/instrumentation/wfc3/performance/throughputs


# In[2]:


madcows1 = fits.open("C:/Users/19133/Documents/cows_official.fits")
madcows = fits.open("C:/Users/19133/Documents/MC2.fits")


# In[3]:


plt.hist(madcows1[1].data['dec'])


# In[4]:


plt.hist(madcows[1].data['dec'])


# In[ ]:


obsids_for_data = [] 
madcow_sources = [] 

failed = [] 

for x in range(len(madcows[1].data['RA'])):
    try: 
        obs_table = Observations.query_object(str(madcows[1].data['RA'][x]) + " " + str(madcows[1].data['DEC'][x]),radius="0.03 deg")
        inds = np.where(obs_table['instrument_name'] == "WFC3/IR")[0]
        if len(inds) > 0:
            for y in obs_table[inds]['obsid']:
                obsids_for_data.append(y)
                madcow_sources.append(x)
        else:
            True 
    except: 
        failed.append(x)


# In[ ]:


obsids_for_data


# In[ ]:


#look at x and continue from there 
#make sure these are in degrees 
#just retry one it failed on? 
#glob images and choose good ones 


# In[ ]:


for x in obsids_for_data:
    data_products = Observations.get_product_list(x)
    manifest = Observations.download_products(data_products, productType="PREVIEW")


# In[16]:


import csv

with open('C:/Users/19133/Documents/MAST_matching_targets_2.csv', 'w', newline='') as csvfile:
    fieldnames = ['obsid','MADCOW_id']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for x in range(len(obsids_for_data)):
        writer.writerow({'obsid':obsids_for_data[x],'MADCOW_id':madcow_sources[x]})


# In[ ]:


#used this in online interface
import csv

with open('C:/Users/19133/Documents/madcows_targets.csv', 'w', newline='') as csvfile:
    fieldnames = ['RA', 'DEC']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for x in range(len(madcows[1].data['RA_SGML'])):
        writer.writerow({'RA': madcows[1].data['RA'][x], 'DEC': madcows[1].data['DEC'][x]})


# In[ ]:


from astropy.table import Table
matches = Table.read("MAST_matching_targets.csv")
obsids = matches['obsid']

for x in [0,1,2]:
    data_products = Observations.get_product_list(str(obsids[x]))
    #print(data_products['dataURI'])
    manifest = Observations.download_products(data_products, productType="SCIENCE")


# In[ ]:


#madcows 2 matching 

#glob jpgs
#go from there, ds9 for best ones 

#then irac stuff 

