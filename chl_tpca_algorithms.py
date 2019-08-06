#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed June 12, 2019

This module contains a functions to calculate the tropical Pacific chlorophyll algorithm (TPCA) concentrations for the satellite sensors; SeaWiFS, MODIS-Aqua and MERIS
Developed for Pittman et al., 2019. JGR: Oceans (2019JC015498)

Notes:
    Numpy has been used over pure python or the math library for more efficient implementations. For increased efficiency, dask (dask.array) can be used instead of numpy.
    These functions are designed for Level 3 Mapped products, which may not precisely reproduce the implementations available on the NASA website (https://oceandata.sci.gsfc.nasa.gov/). 

General functions include: 
    blended_chl
    calculate_chl_ocx
    calculate_chl_ci
    
Sensor specific functions include:
    calculate_seawifs_chl
    calcuate_modis_chl
    calculate_meris_chl

@author: Nicholas.Pittman
@email: Nic.Pittman@utas.edu.au
@position: PhD Candidate; Biogeochemistry, remote sensing, oceanography, Tropical Pacific
@affiliation1: Institute of Marine and Antarctic Studies, University of Tasmania
@affiliation2: Australian Research Council Centre of Excellence for Climate Extremes.

References:
    Tropical Pacific Chlorophyll Algorithm (2019JC015498)
    'An assessment and improvement of satellite ocean color algorithms for the tropical Pacific Ocean'
    Pittman, N., Strutton, P., Matear, R., Johnson, R., (2019)
    Submitted: Journal of Geophysical Research: Oceans
    
    CI and blending window
    Hu, C., Lee, Z., and Franz, B. (2012). 
    Chlorophyll a algorithms for oligotrophic oceans: A novel approach based on three-band reflectance difference.
    Journal of Geophysical Research: Oceans 117.
    
    Traditional OCx Algorithm
    O’Reilly, J.E., Maritorena, S., Mitchell, B.G., Siegel, D.A., 
    Carder, K.L., Garver, S.A., Kahru, M., and McClain, C. (1998).
    Ocean color chlorophyll algorithms for SeaWiFS. 
    Journal of Geophysical Research: Oceans 103, 24937–24953.
"""

import numpy as np       #Version: '1.16.1'
#import dask.array as np #Version: '1.0.0'

def blended_chl(chl_ci,chl_ocx,l=0.15,h=0.2):  #Default blending window of 0.15 to 0.2
    """A general Chl algorithm blending function between Chl_CI to Chl_OCx"""
    upper=h #0.2
    lower=l #0.15
    alpha=(chl_ci-lower)/(upper-lower)
    beta=(upper-chl_ci)/(upper-lower)
    chl = (alpha*chl_ocx)+beta*chl_ci
    chl = np.where(chl_ci < lower, chl_ci, chl)
    chl = np.where(chl_ci > upper, chl_ocx, chl)
    return chl

def calculate_chl_ocx(ocx_poly,lmbr):
    """A general Chl_OCx algorithm for fourth order polynomials (O'Reilly et al., 1998)"""
    chl_ocx=10**(ocx_poly[0]+(ocx_poly[1]*lmbr**1)+(ocx_poly[2]*lmbr**2)+(ocx_poly[3]*lmbr**3)+(ocx_poly[4]*lmbr**4))
    return chl_ocx
    
def calculate_chl_ci(ci_poly,CI):
    """A general Chl_CI algorithm for a linear polynomial (Hu et al., 2012)"""
    chl_ci=10**(ci_poly[0]+ci_poly[1]*CI)
    return chl_ci

def calculate_seawifs_chl(r443,r490,r510,r555,r670,ocx_poly=[0.3255,-2.7677,2.4409,-1.1288,-0.4990],ci_poly=[-0.4909, 191.6590],l=0,h=0.5):
    """
    Given:
        SeaWiFS RRS values for 443,490,510,555,670
        OCx Polynomial
        CI Polynomial
        l - Low blending cutoff
        h - High blending cutoff 
        
        
    Calculate:
        Calculate Chl OCx
        Calculate Chl CI
        Return blended chlorophyll esimate (Default product is the Pittman et al., 2019 TPCA)
    
    Usage:
        #Calculate TPCA SeaWiFS Chlorophyll
        calculate_seawifs_chl(rrs443,rrs490,rrs510,rrs555,rrs670, poly_ocx=[0.3255,-2.7677,2.4409,-1.1288,-0.4990],l=0,h=0.5)
        
        #Calculate NASA SeaWiFS Chlorophyll
        calculate_seawifs_chl(rrs443,rrs490,rrs510,rrs555,rrs670, poly_ocx=[0.3272,-2.9940, 2.7218,-1.2259,-0.5683],l=0.15,h=0.2)
  
    Polynomials:
        TPCA SeaWiFS: [0.3255,-2.7677,2.4409,-1.1288,-0.4990], l=0,h=0.5 (OCx)
        NASA SeaWiFS: [0.3272,-2.9940, 2.7218,-1.2259,-0.5683], l=0.15,h=0.2 (OCx)
        Hu2012: [-0.4909, 191.6590] #Same for both algorithms (CI)
    """
    
    #Calculate Chl OCX (O'Reilly et al., 1998)
    mbr=np.maximum(r443/r555,r490/r555)
    mbr=np.maximum(mbr,r510/r555) #Calculate max band ratio

    lmbr=np.log10(mbr)
    chl_ocx=calculate_chl_ocx(ocx_poly,lmbr)
        
    #Calculate Chl CI (Hu et al., 2012)
    CI=r555-(r443+(555-443)/(670-443)*(r670-r443))
    chl_ci=calculate_chl_ci(ci_poly,CI)
    
    #Blending between Chl_CI to Chl_OCx
    blended=blended_chl(chl_ci,chl_ocx,h=h,l=l) 
    return blended


def calculate_modis_chl(r443,r488,r547,r667,ocx_poly=[0.3272,-2.9940,2.7218,-1.2259,-0.5683],ci_poly=[-0.4909, 191.6590],l=0,h=0.2):
    """
    Given:
        MODIS-Aqua RRS values for 443,488,547,667
        OCx Polynomial
        CI Polynomial
        l - Low blending cutoff
        h - High blending cutoff 

        
    Calculate:
        Calculate Chl OCx
        Calculate Chl CI
        Return blended chlorophyll esimate (Default product is the Pittman et al., 2019 TPCA)
    
    Usage:
        #Calculate TPCA MODIS-Aqua Chlorophyll
        calc_modis_chl(rrs443,rrs490,rrs510,rrs555,rrs670, poly_ocx=[0.3272,-2.9940, 2.7218,-1.2259,-0.5683],l=0,h=0.2)
        
        #Calculate NASA MODIS-Aqua  Chlorophyll    
        calculate_modis_chl(rrs443,rrs490,rrs510,rrs555,rrs670, poly_ocx=[0.2424,-2.7423,1.8017,0.0015,-1.2280],l=0.15,h=0.2)

    Polynomials:
        TPCA MODIS-Aqua: [0.3272,-2.9940, 2.7218,-1.2259,-0.5683], l=0,h=0.5 (OCx)
        NASA MODIS-Aqua: [0.2424,-2.7423,1.8017,0.0015,-1.2280], l=0.15,h=0.2 (OCx)
        Hu2012: [-0.4909, 191.6590] #Same for both algorithms (CI)
    """
        
    #Calculate Chl OCX (O'Reilly et al., 1998)
    mbr=np.maximum(r443/r547,r488/r547) #Calculate max band ratio
    lmbr=np.log10(mbr)
    chl_ocx=calculate_chl_ocx(ocx_poly,lmbr)
    #Calculate Chl CI (Hu et al., 2012)
    CI=r547-(r443+(547-443)/(667-443)*(r667-r443))
    chl_ci=calculate_chl_ci(ci_poly,CI)
    
    blended=blended_chl(chl_ci,chl_ocx,h=h,l=l) #Blending cutoff
    return blended


def calculate_meris_chl(r443,r490,r510,r560,r665,ocx_poly=[0.3255,-2.7677, 2.4409,-1.1288,-0.4990],ci_poly=[-0.4909, 191.6590],l=0.15,h=0.2):
    """
    Given:
        MERIS RRS values for 443,490,510,560,665
        OCx Polynomial
        CI Polynomial
        l - Low blending cutoff
        h - High blending cutoff 
        
    Calculate:
        Calculate Chl OCx
        Calculate Chl CI
        Return blended chlorophyll esimate (Default product is the Pittman et al., 2019 TPCA)
        MERIS TPCA is the default NASA implementation.
    
    Usage:
        #Calculate NASA / TPCA MODIS-Aqua  Chlorophyll (identical)
        calculate_meris_chl(rrs443,rrs490,rrs510,rrs560,rrs665,poly_ocx=[0.3255,-2.7677, 2.4409,-1.1288,-0.4990],l=0.15,h=0.2)

    Polynomials:
        NASA MERIS: [0.3255,-2.7677, 2.4409,-1.1288,-0.4990], l=0.15,h=0.2 (OCx)
        Hu2012: [-0.4909, 191.6590] #Same for both algorithms (CI)
    """
    
    #Calculate Chl OCX (O'Reilly et al., 1998)
    mbr=np.maximum(r443/r560,r490/r560)
    mbr=np.maximum(mbr,r510/r560) #Calculate max band ratio
    lmbr=np.log10(mbr)
    chl_ocx=calculate_chl_ocx(ocx_poly,lmbr)
    
    #Calculate Chl CI (Hu et al., 2012)
    CI=r560-(r443+(560-443)/(665-443)*(r665-r443))
    chl_ci=calculate_chl_ci(ci_poly,CI)
    
    blended=blended_chl(chl_ci,chl_ocx,h=h,l=l) #Blending cutoff
    return blended


if __name__ == '__main__':
    pass
    
