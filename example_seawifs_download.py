#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A simple test script to download the 2001-01-01 SeaWiFS RRS data, process with the chl_tpca_algorithms functions, plot and compare to the NASA chlor_a data.

@author: Nicholas.Pittman
@email: Nic.Pittman@utas.edu.au
@position: PhD Candidate; Biogeochemistry, remote sensing, oceanography, Tropical Pacific
@affiliation1: Institute of Marine and Antarctic Studies, University of Tasmania
@affiliation2: Australian Research Council Centre of Excellence for Climate Extremes.
"""

from chl_tpca_algorithms import calculate_seawifs_chl

import requests                      #Version '2.19.1'
import xarray as xr                  #Version '0.11.3'
import numpy as np                   #Version '1.16.1' (Not used here, but in chl_tpca_algorithms)
import matplotlib.pyplot as plt      #Version '3.0.0'
import os

##################################################### General functions, Download and cutout
def downloader(urls):
    """" Function designed to download satellite data from https://oceandata.sci.gsfc.nasa.gov/"""
    
    def exists(fileloc):
        try:
            xr.open_dataset(fileloc)
            return True
        except:
            return False
    
    path='seawifs_data'
    if not os.path.isdir(path):
        print('Creating directory: ',path)
        os.makedirs(path)    
        
    file_locations=[]
    for url in urls:
        while True:
            #Download the files to their file name in the directory we just created.
            #ie: seawifs_data/S2000001.L3m_DAY_RRS_Rrs_443_9km.nc
            fileloc=path+'/'+url.split('/')[-1]
            if exists(fileloc):
                print('Exists: ',fileloc)
                file_locations.append(fileloc)
                break
            r = requests.get(url)#,timeout=20)
            with open(fileloc, 'wb') as f:
                f.write(r.content)

            #Ensure that the file actually downloaded, this can fail sometimes for some reason, maybe too soon.
            #time.sleep(1) #Can fail sometimes so maybe this is a fix
            if (r.status_code==200) & (exists(fileloc)==True):
                print('Downloaded:',fileloc)
                file_locations.append(fileloc)
                break
            else:
                print('Download failed:', fileloc,'status:',r.status_code)
                
    return file_locations

def cut_tropical_pacific(chl_dataset):
    """Cut the tropical pacific out of the global dataset, move the dateline 180° and clean up a little bit"""
    
    chl_dataset= chl_dataset.assign_coords(lon=(chl_dataset.lon % 360)).roll(lon=(chl_dataset.dims['lon'] // 2),roll_coords=True) #Change the dateline to the otherside of the world
    chl_trop_pac=chl_dataset.sel(lat=slice(20,-20),lon=slice(120,290))
    #print('Subset Size: ',chl_trop_pac.nbytes / 1e9, 'GB')
    return chl_trop_pac

##################################################### Download, process and run 

urls=['https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S2000001.L3m_DAY_RRS_Rrs_443_9km.nc',
      'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S2000001.L3m_DAY_RRS_Rrs_490_9km.nc',
      'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S2000001.L3m_DAY_RRS_Rrs_510_9km.nc',
      'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S2000001.L3m_DAY_RRS_Rrs_555_9km.nc',
      'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S2000001.L3m_DAY_RRS_Rrs_670_9km.nc',
      'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S2000001.L3m_DAY_CHL_chlor_a_9km.nc']


#Download files in a consistent way    
file_locations=downloader(urls)

#Open the files with xarray, and use the cutout function to get the tropical Pacific out.
rrs = cut_tropical_pacific(xr.open_mfdataset(file_locations[0:5]))
nasa_chl_a = cut_tropical_pacific(xr.open_dataset(file_locations[5])).chlor_a

#Process the files
chl_tpca=xr.apply_ufunc(calculate_seawifs_chl, rrs.Rrs_443,rrs.Rrs_490,rrs.Rrs_510,rrs.Rrs_555,rrs.Rrs_670, dask='allowed')

##################################################### Plot some maps

#Plot TPCA map
pallette=plt.cm.viridis
pallette.set_bad(color='gray') #Just so we can have a gray rather than white background.
plt.figure(figsize=(10,5))
chl_tpca.plot(vmin=0,vmax=0.5,cmap=pallette)
plt.title('Tropical Pacific Chlorophyll Algorithm (mg m$^{-3}$, SeaWiFS: 2000-01-01)')
plt.ylabel('Latitude')
plt.xlabel('Longitude (inverse)')
plt.show()

#Plot NASA Chlor_a map
plt.figure(figsize=(10,5))
nasa_chl_a.plot(vmin=0,vmax=0.5,cmap=pallette)
plt.title('NASA Chlor a (mg m$^{-3}$, SeaWiFS: 2000-01-01)')
plt.ylabel('Latitude')
plt.xlabel('Longitude (inverse)')
plt.show()

#Plot Differences between TPCA and NASA
pallette_bwr=plt.cm.bwr
pallette_bwr.set_bad(color='gray') #Just so we can have a gray rather than white background.
plt.figure(figsize=(10,5))
difference=nasa_chl_a-chl_tpca
difference.plot(vmin=-0.05,vmax=0.05,cmap=pallette_bwr)
plt.title('NASA Chlor a - TPCA (mg m$^{-3}$, SeaWiFS: 2000-01-01)')
plt.ylabel('Latitude')
plt.xlabel('Longitude (inverse)')
plt.show()
