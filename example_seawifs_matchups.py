#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A short script to use the SeaWiFS matchup data to reproduce the TPCA and reproduce Pittman et al. 2019 statistics and Figure 8a.

Note***
The diagnostics produced by this script do not identically reproduce Table 3, SeaWiFS rank 2.
The Rrs values provided in the matchup databases are an average of the 45 pixel matchup window. 
Table 3 was produced by calculating chl for each of the 45 pixels, and then averaging the 45 chlorophyll concentrations.
This produces slightly different values than seen in the paper. 

@author: Nicholas.Pittman
@email: Nic.Pittman@utas.edu.au
@position: PhD Candidate; Biogeochemistry, remote sensing, oceanography, Tropical Pacific
@affiliation1: Institute of Marine and Antarctic Studies, University of Tasmania
@affiliation2: Australian Research Council Centre of Excellence for Climate Extremes.
"""

from chl_tpca_algorithms import calculate_seawifs_chl
import pandas as pd
from chl_statistics import plot_linear_trend, check_bias

seawifs_matchups_path='tropical_pacific_matchups/seawifs_matchups.csv'
seawifs_matchups=pd.read_csv(seawifs_matchups_path)

seawifs_matchups=seawifs_matchups[seawifs_matchups.validation_set==False]
print('Columns in the matchup files:',seawifs_matchups.columns.values) #So we can see what is in the data file. 

#Calculating the TPCA is this easy!
tpca_chl=calculate_seawifs_chl(seawifs_matchups.rrs443,seawifs_matchups.rrs490,seawifs_matchups.rrs510,seawifs_matchups.rrs555,seawifs_matchups.rrs670)
tpca_chl_r=calculate_seawifs_chl(seawifs_matchups.rrs443,seawifs_matchups.rrs490,seawifs_matchups.rrs510,seawifs_matchups.rrs555,seawifs_matchups.rrs670,ocx_poly=[0.3272,-2.9940, 2.7218,-1.2259,-0.5683], l=0.15,h=0.2)


#A check which can reproduce some statistics and Figure 8a in Pittman et al., 2019. 
check_bias(seawifs_matchups.in_situ_chl, tpca_chl,seawifs_matchups.NASA_chlor_a)
