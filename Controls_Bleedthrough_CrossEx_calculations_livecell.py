# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 14:28:10 2021

@author: guada
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# %% general parameters for all figures
mpl.style.use('default')
mpl.rcParams['axes.linewidth'] = 7  # set the value globally
mpl.rcParams['xtick.major.size'] = 20
mpl.rcParams['xtick.major.width'] = 7
mpl.rcParams['xtick.minor.size'] = 10
mpl.rcParams['xtick.minor.width'] = 7
mpl.rcParams['ytick.major.size'] = 20
mpl.rcParams['ytick.major.width'] = 7
mpl.rcParams['ytick.labelsize'] = 50
mpl.rcParams['xtick.labelsize'] = 50
mpl.rcParams['ytick.minor.size'] = 10
mpl.rcParams['ytick.minor.width'] = 7
mpl.rcParams['font.size'] = 55
mpl.rcParams['font.sans-serif'] = 'Arial'
#%% data import 
df=pd.read_csv("FPcontrols_livecell.csv")
colors=pd.read_csv("GSs_color_scheme.csv")
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
    fig1,ax1=plt.subplots(1,1,figsize=(10,10))
    # sns.regplot(x,y,data=selected)
    plt.scatter(x,y,data=selected)
    plt.plot(selected[x],np.polyval(p,selected[x]),color='red')
    ax1.set_ylabel(y,size=60)
    ax1.set_xlabel(x,size=60)
    # plt.ylim(2500,8000)
    plt.title('y={:.2f}x + {:.2f} \nslope \u03C3 \u00B1 {:.3f} \n  b-int \u03C3 \u00B1 {:.0f}'.format(slope, intercept,slope_std_err,int_std_err),horizontalalignment='center',x=0.3, y=0.7,size=35)
    plt.locator_params(axis='x',nbins=4)
    plt.locator_params(axis='y',nbins=4)
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    if construct=='mTQ':
        plt.suptitle('Bleedthrough', size=60)
    else:
        plt.suptitle('Cross excitation', size=60)

#%% fit for getting bleed through

## IMPORTANT notes
## for cells expressing only mTurquoise2 [select data that contains only mtq2]
## plotting donor emission (under FRET) on the x-axis (ch1)
## plotting acceptor emission under FRET on the y-axis (ch2)

plot_FP(df,'mTQ',"D",'A')
plt.ylabel('acceptor emission \n(ch2)',size=60)
plt.xlabel("donor emission \n(ch1)",size=60)
plt.savefig('bleedthrough.svg', format="svg",bbox_inches='tight')
#%% fit for getting cross excitation

## IMPORTANT notes
## for cells expressing only mNeonGreen [select data that contains only mNG]
## plotting direct acceptor on the x-axis (ch3)
## plotting acceptor emission under FRET on the y-axis (ch2)
# plot_FP(df,'mNG',"ch3",'ch2')
plot_FP(df,'mNG',"directA",'A')
plt.ylabel('acceptor emission \n(ch2)',size=60)
plt.xlabel("direct acceptor emission \n(ch3)",size=60)
plt.savefig('crossexcitation.svg', format="svg",bbox_inches='tight')
#%% function to calculate slopes (bleed through and cross excitation) 
def find_slope(df,construct,x,y):
    selected=df[df['construct']==construct]
    p,cov=np.polyfit(selected[x], selected[y],1,cov=True)
    slope=p[0]
    intercept=p[1]
    std_err=np.sqrt(np.diag(cov))
    slope_std_err= std_err[0]
    int_std_err = std_err[1]
    return(slope)
#%% applying correction factors
bleedthrough= 0.39 #find_slope(df,'mTQ',"D",'A')
crossexcitation= 0.1 #find_slope(df,'mNG',"directA",'A')
df['A_corrected']=df['A']-(df['D']*bleedthrough)-(df['directA']*crossexcitation)
df['Ef']=df['A']/(df['D']+df['A_corrected'])
#%% defining plot_GSs function: violins of all the GS linkers 
def plot_GSs(df1,Ef_col):
    allMeans = pd.DataFrame()
    GSs = ['GS0','GS16','GS32','GS48']
    fig,ax = plt.subplots(1,1, figsize=[10,10], sharex=True, sharey=True)
    means=np.array([])
    errs=np.array([])
    ax.grid()
    for protIdx,prot in enumerate(GSs):
        sliced = df1[df1.construct==prot]
        color = colors[colors['construct']==prot]['color']
        N_res = int(prot[2:])*2
        parts = ax.violinplot(sliced[Ef_col],positions=[N_res],
                              widths=20,showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(color)
            pc.set_edgecolor('black')
            pc.set_lw(3)
            pc.set_alpha(1)
        bot, quartile1, medians, quartile3, top = np.percentile(sliced[Ef_col], [5, 25, 50, 75, 95])
        allMeans = allMeans.append({'color':color.values[0],
                                    'prot':prot,'q1':quartile1,'median':medians,'q3':quartile3},ignore_index=True)
        ax.vlines(N_res, quartile1, quartile3, color='r', linestyle='-', lw=10)
        ax.vlines(N_res, bot, top, color='r', linestyle='-', lw=3)
        ax.scatter(N_res,medians,alpha=1,marker='s',
                   c='w',s=80,zorder=3,linewidth=2,edgecolor='k')
        means = np.append(means,medians)
        errs = np.append(errs,quartile3-quartile1)
    Ef_GS_x=np.array([int(x[2:])*2 for x in iter(GSs)])
    fit,cov = np.polyfit(Ef_GS_x,means,1,cov=True,w=1/errs)
    fit_err = np.sqrt(np.diag(cov))
    Ef_GS_x=np.append(Ef_GS_x,200)
    Ef_GS_y=Ef_GS_x*fit[0]+fit[1]
    Ef_GS_y_top = Ef_GS_x*(fit[0]+fit_err[0])+fit[1]+fit_err[1]
    Ef_GS_y_bot = Ef_GS_x*(fit[0]-fit_err[0])+fit[1]-fit_err[1]
    ax.plot(Ef_GS_x,Ef_GS_y,'--',c='cadetblue',zorder=3,lw=5)
    ax.fill_between(Ef_GS_x,Ef_GS_y_bot,Ef_GS_y_top,color='cadetblue',alpha=0.2,zorder=3)
    ax.set_xlabel('$N_{residues}$',fontsize=60)
    ax.set_ylabel('$E_f^{cell}$',fontsize=60)
    ax.set_ylim(0,1)
    ax.set_xlim(-10,170)
    ax.set_xticks([0,50,100,150])
    plt.savefig('GSs_trends.svg', format="svg",bbox_inches='tight', dpi=1200)

#%% plotting GS series
GSs_df=df[df['construct'].str.startswith('GS')]

plot_GSs(GSs_df,'Ef_dirty')
plot_GSs(GSs_df,'Ef')


