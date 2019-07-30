#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Includes:
    plot_linear_trend
        A linear trend plotter which returns linear statistics.
    check_bias
        Bias diagnostics including: % Wins, mean and median log bias+absolute error.
Small function to plot linear trends of data (Specifically low chlorophyll concentrations)

@author: npittman
"""

import numpy as np                   #Version '1.16.1'
import matplotlib.pyplot as plt      #Version '3.0.0'
from scipy.stats import linregress   #Version '1.1.0'


def plot_linear_trend(x,y,title,xlab,ylab,diagonal_line=1,printer=1,plot=1,logspace=0,color=np.ndarray(0),zx=np.ndarray(0),zy=np.ndarray(0),trendline=1):
    """"
    Designed for chlorophyll trends <1mg/m3 
    
    Given an x and y (same size)
    + Title, X Label, Y Label
    
    if Z, will plot 
    Optional: Diagonal_line=0, default on.
    Optional: Printer=0, default on
    
    Returns linear statistics. 
    """
    
    x=np.ravel(x)
    y=np.ravel(y)
    mask=~np.isnan(x)
    x=x[mask]
    y=y[mask]
    
    slope, intercept, r_value, p_value, std_err = linregress(x,y)
    if plot==0:
         return slope, intercept, r_value, p_value, std_err 
    
    mn=0
    mx=1
    x1=np.linspace(mn,mx,300)
    y1=slope*x1+intercept
    plt.figure(figsize=(10,10))
    if zx.size==0:
        if color.size==0:
            plt.scatter(x,y,c='k'),plt.xlim([0,0.5]),plt.ylim([0,0.5]) #old blue #061264
        else:
            plt.scatter(x,y,c=color,cmap='coolwarm',vmin=0,vmax=0.5),plt.colorbar(),plt.xlim([0,0.6]),plt.ylim([0,0.6])
    else:
        plt.scatter(zx,zy,c='none',marker='o',edgecolors='k')
        plt.scatter(x,y,c='k'),plt.xlim([0,0.5]),plt.ylim([0,0.5])
        
        slopea, intercepta, r_valuea, p_valuea, std_erra = linregress(zx,zy)
        x2=np.linspace(mn,mx,300)
        y2=slopea*x2+intercepta
        if trendline==1:
            plt.plot(x2,y2,':k',linewidth=2)  
        
    if logspace==1:
        plt.gca().set_xscale('log'),plt.gca().set_yscale('log')
        plt.xlim([0,1]),plt.ylim([0,1])
    
    plt.title(title,fontsize=32),plt.xlabel(xlab,fontsize=28),plt.ylabel(ylab,fontsize=28)
    plt.tick_params(axis='both', which='major', labelsize=24)    #fontsize for ticks

    if diagonal_line==1:
        plt.plot(np.arange(0,1,0.001),np.arange(0,1,0.001),'--r')
    if trendline==1:
        plt.plot(x1,y1,'-g',linewidth=3)  
    #plt.savefig('../results/figures/'+title+'.png', format='png', dpi=300,transparent=True)   
    plt.show()
    if printer==1:
        print(title,'Slope: ',slope,' Intercept: ',intercept,' R2: ',r_value)
     
    return slope, intercept, r_value, p_value, std_err



def check_bias(in_situ,model_a,model_b,plot=1):
    a_bias=((model_a-in_situ)/in_situ)*100
    a_abs_error=abs(a_bias)
    
    b_bias=((model_b-in_situ)/in_situ)*100
    b_abs_error=abs(b_bias)
    
    a_wins,b_wins,draw=0,0,0
    try:
        for i in range(len(in_situ)):
            if a_abs_error.iloc[i]<b_abs_error.iloc[i]:
                a_wins+=1
            elif a_abs_error.iloc[i]>b_abs_error.iloc[i]:
                b_wins+=1
            else:
                draw+=1
    except:
        for i in range(len(in_situ)):
            if a_abs_error[i]<b_abs_error[i]:
                a_wins+=1
            elif a_abs_error[i]>b_abs_error[i]:
                b_wins+=1
            else:
                draw+=1
    total=a_wins+b_wins        
    #print(a_wins,b_wins,draw) 
    try:
        a_win_percent=(a_wins/total)*100
    except:
        print('nothing?')
        
    try:
        b_win_percent=(b_wins/total)*100
    except:
        print('nothing?')
        
    new_bias=10**(np.nanmean(np.log10(model_a)-np.log10(in_situ)))
    nasa_bias=10**(np.nanmean(np.log10(model_b)-np.log10(in_situ)))
    
    new_med_bias=10**(np.nanmedian(np.log10(model_a)-np.log10(in_situ)))
    nasa_med_bias=10**(np.nanmedian(np.log10(model_b)-np.log10(in_situ)))
    
    
    new_mae=10**(np.nanmean(abs(np.log10(model_a)-np.log10(in_situ))))
    nasa_mae=10**(np.nanmean(abs(np.log10(model_b)-np.log10(in_situ))))
        
    new_med_ae=10**(np.nanmedian(abs(np.log10(model_a)-np.log10(in_situ))))
    nasa_med_ae=10**(np.nanmedian(abs(np.log10(model_b)-np.log10(in_situ))))
    
    slope_a,int_a,r2_a,_,_=plot_linear_trend(in_situ,model_a,'New Algorithm vs Observed','In situ $mg/m^3$','New Algorithm Chlor a',plot=0)
    slope_b,int_b,r2_b,_,_=plot_linear_trend(in_situ,model_b,'NASA Algorithm vs Observed','In situ $mg/m^3$','NASA Chlor a',plot=0)
    
    plot_linear_trend(in_situ,model_a,'Both Algorithms vs Observed','In situ $mg/m^3$','New Algorithm Chlor a',plot=plot,zx=in_situ,zy=model_b)
    
    if plot==1:
        print('Model A')   
        print('Model A was better:',np.round(a_win_percent,3),'% of observations')
        #print('Model A Mean Log Bias: ',np.round(new_bias,3))
        print('Model A Median Log Bias: ',np.round(new_med_bias,3))
        #print('Model A Mean AE: ',np.round(new_mae,3))
        print('Model A Median AE: ',np.round(new_med_ae,3))
        print('Model A Slope: ',np.round(slope_a,3))
        print('Model A Intercept: ',np.round(int_a,3))
        print('Model A R2: ', np.round(r2_a,3))

        print('')
        print('Model B')
        print('Model B was better:',np.round(b_win_percent,3),'% of observations')
        #print('Model B Mean Log Bias: ',np.round(nasa_bias,3))
        print('Model B Median Log Bias: ',np.round(nasa_med_bias,3))
        #print('Model B Mean AE: ',np.round(nasa_mae,3))
        print('Model B Median AE: ',np.round(nasa_med_ae,3))
        print('Model B Slope: ',np.round(slope_b,3))
        print('Model B Intercept: ',np.round(int_b,3))
        print('Model B R2: ', np.round(r2_b,3))
    
    