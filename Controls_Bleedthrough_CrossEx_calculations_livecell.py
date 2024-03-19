# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 14:28:10 2021

@author: guada
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#%% data import 
df=pd.read_csv("FPcontrols_livecell.csv")
#%% merging to construct and condition
conditions=pd.read_csv("FPcontrols_plate_layout.csv")
df=pd.merge(df,conditions,on=['construct','well'],how='outer').dropna()
#%% defining plot_FP function: will slice data, do linear fits between two channels
def plot_FP(df,construct,x,y):
    selected=df[df['construct']==construct]
    p,cov=np.polyfit(selected[x], selected[y],1,cov=True)
    slope=p[0]
    intercept=p[1]
    std_err=np.sqrt(np.diag(cov))
    slope_std_err= std_err[0]
    int_std_err = std_err[1]
    fig1,ax1=plt.subplots(1,1,figsize=(4,4))
    # sns.regplot(x,y,data=selected)
    plt.scatter(x,y,data=selected)
    plt.plot(selected[x],np.polyval(p,selected[x]),color='red')
    ax1.set_ylabel(y,size=18)
    ax1.set_xlabel(x,size=18)
    # plt.ylim(2500,8000)
    plt.title('y={:.2f}x + {:.2f} \nslope \u03C3 \u00B1 {:.3f} \n  b-int \u03C3 \u00B1 {:.0f}'.format(slope, intercept,slope_std_err,int_std_err),horizontalalignment='center',x=0.3, y=0.7)
    plt.locator_params(axis='x',nbins=4)
    plt.locator_params(axis='y',nbins=4)
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    if construct=='mTQ2':
        plt.suptitle('Bleedthrough', size=20)
    else:
        plt.suptitle('Cross excitation', size=20)

#%% fit for getting bleed through

## IMPORTANT notes
## for cells expressing only mTurquoise2 [select data that contains only mtq2]
## plotting donor emission (under FRET) on the x-axis (ch1)
## plotting acceptor emission under FRET on the y-axis (ch2)

plot_FP(df,'mTQ2',"ch1",'ch2')
plt.ylabel('acceptor emission \n(ch2)',size=18)
plt.xlabel("donor emission \n(ch1)",size=18)
plt.savefig('bleedthrough.svg', format="svg",bbox_inches='tight')
#%% fit for getting cross excitation

## IMPORTANT notes
## for cells expressing only mNeonGreen [select data that contains only mNG]
## plotting direct acceptor on the x-axis (ch3)
## plotting acceptor emission under FRET on the y-axis (ch2)

plot_FP(df,'mNG',"ch3",'ch2')
plt.ylabel('acceptor emission \n(ch2)',size=18)
plt.xlabel("direct acceptor emission \n(ch3)",size=18)
plt.savefig('crossexcitation.svg', format="svg",bbox_inches='tight')
