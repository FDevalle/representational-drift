#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 16:34:35 2022

@author: federico
"""


import numpy as np
import matplotlib.pyplot as plt

import set_figure
from matplotlib.pyplot import cm
def plotting_statistics(dataexp,datamodel,network_data,model_avg_act_vec_spat,active_mice,figure_name,ndays=8,plotting_exp=True,
                                    xsize=1.6,ysize=3,legend=False):
    """
    

    Parameters
    ----------
    dataexp : array
        Experimental data of dimension (3*).
    datamodel : TYPE
        DESCRIPTION.
    network_data : TYPE
        DESCRIPTION.
    model_avg_act_vec_spat : TYPE
        DESCRIPTION.
    active_mice : TYPE
        DESCRIPTION.
    figure_name : TYPE
        DESCRIPTION.
    ndays : TYPE, optional
        DESCRIPTION. The default is 8.
    plotting_exp : TYPE, optional
        DESCRIPTION. The default is True.
    xsize : TYPE, optional
        DESCRIPTION. The default is 1.6.
    ysize : TYPE, optional
        DESCRIPTION. The default is 3.
    legend : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    data_total=dataexp  
    avg_solutions_spat=datamodel
    color_data='darkblue'
    color_network='limegreen'
    color_network_sym='coral'
    color_stat='red'
    color_stat_spat='red'
    
    scatt_net_col='limegreen'
    scatt_net_col_sym='coral'
    
    scatt_data_col='darkblue'
    stat_width=1
    w_line=1
    ms=3
    ms2=5
    xlab=[1,8]
    xlabdelta=[0,2,4,6]
    w_line_sh=1
    alpha_sh=0.8
    padx=-1.5
    st_sh='--'
    #active_mice=dataexp[7:15,:]
    #color = cm.Reds(np.linspace(0.2, 1, avg_solutions_spat.shape[1]))

    fig,ax=set_figure.set_figure(3,1,xsize,ysize)
    
    
    #fig.suptitle(r"$p=%1.2f$, $\tilde{\sigma}=%1.2f$, $\tilde{\theta}=%1.2f$, $p_s=%1.2f$" 
    #             %(popt0[1],popt0[2],popt0[0],popt0[4]),fontsize=16,y=0.9)
    
    #ax[0,0].set_title("Active cells",y=1.3)
    
    #ndays=8
    sessions=np.arange(1,ndays+1)
    session2=np.arange(1,ndays) 
    local_overl=np.concatenate(([1],avg_solutions_spat[:ndays-1]))

    ax[1].plot(session2,avg_solutions_spat[3*ndays:],color=color_stat,label='stat. model',zorder=3,linewidth=stat_width)
    ax[2].plot(np.arange(1,ndays),local_overl[1:],color=color_stat,label='stat. model',zorder=3,linewidth=stat_width)
    ax[0].plot(np.arange(1,ndays+1),avg_solutions_spat[2*ndays-1:3*ndays-1],color=color_stat,label='stat. model',zorder=3,linewidth=stat_width)


    #data_time_full=[1.,         0.70779221, 0.57272727, 0.50194805, 0.45194805, 0.4025974,
    # 0.37077922, 0.35064935]
    if plotting_exp:
        ax[1].errorbar(session2,data_total[3*ndays:,0],yerr=data_total[3*ndays,1],fmt='-o',color=color_data,markersize=ms,mfc=scatt_data_col,linewidth=w_line,label=' exp. data')
    ax[1].errorbar(session2,network_data[ndays+1:2*ndays,0],yerr=network_data[ndays+1:2*ndays,1],fmt='-*',color=color_network,markersize=ms2,mfc=scatt_net_col,linewidth=w_line,label=' network model')
    
    #ax[1].plot(session2,avg_solutions_spat[1:8],color=color_stat_spat,label='stat. model',zorder=3,linewidth=stat_width)
    
    
    #ax[1].fill_between(session2,data_95m_spat[1:8],data_95p_spat[1:8],facecolor=color_stat_spat,alpha=0.15)
    #ax[1,0].plot(session2,avg_surv2[1:])
    #ax[1,0].plot(session2,avg_surv3[1:],color='orange',linewidth=1)
    ax[1].set_ylim([0,1.05])
    ax[1].set_yticks([0,1])
    ax[1].set_yticklabels([0,1])
    
    ax[1].set_ylabel('Survival')
    ax[1].set_xticks([1,7])
    
    ax[1].set_xticklabels([1,7])
    #ax[1].legend(frameon=False,fontsize=6)
    #ax[1,0].set_xlabel(r'$\Delta t$ (sessions)',labelpad=padx)
    
    axins = ax[0].inset_axes([0.6,0.7,0.4,0.4])
    axins.errorbar(sessions,active_mice[:,0],yerr=active_mice[:,1],fmt='-o',color=color_data,markersize=2,mfc=scatt_data_col,linewidth=w_line)
    axins.errorbar(sessions,network_data[-8:,0],yerr=network_data[-8:,1],fmt='-*',color=color_network,markersize=3,mfc=scatt_net_col,linewidth=w_line)
    # #axins.plot(sessions,model_avg_act_vec,'-',color=color_stat)
    axins.plot(sessions,model_avg_act_vec_spat,'-',color=color_stat_spat,zorder=3,linewidth=stat_width)
    axins.spines['right'].set_visible(False)
    axins.spines['top'].set_visible(False)
    
    
    axins.set_ylim([0.,0.7])
    axins.set_yticks([0.,0.6])
    
    axins.set_yticklabels([0,60],fontsize=6)
    daysnum=[1,8]
    axins.set_xticks(daysnum)
    axins.set_xticklabels(daysnum,fontsize=6)
    axins.set_xlabel('Session',labelpad=-5,fontsize=6)
    axins.set_ylabel("Active (\%)",fontsize=6,labelpad=0)
    axins.tick_params('x',pad=1,length=2)
    axins.tick_params('y',pad=1,length=2)
    
    
    
    ##third plot
    #print(average_hamming.shape)
    data_mice=np.row_stack(([1,0],data_total[:ndays-1,:]))
    if plotting_exp:

        ax[2].errorbar(np.arange(1,ndays),data_mice[1:,0],yerr=data_mice[1:,1],fmt='-o',color=color_data,markersize=ms,mfc=scatt_data_col,linewidth=w_line,label=' exp. data')
    #ax[2].plot(np.arange(1,ndays),data_specific[ndays:2*ndays-1,:])
    ax[2].errorbar(np.arange(1,ndays),network_data[2*ndays:3*ndays-1,0],yerr=network_data[2*ndays:3*ndays-1,1],fmt='-*',color=color_network,markersize=ms2,mfc=scatt_net_col,linewidth=w_line,label='network model')
    
    #local_overl=np.concatenate(([1],avg_solutions_spat[8:15]))
    #locm=np.concatenate(([1],data_95m_spat[8:15]))
    #locp=np.concatenate(([1],data_95p_spat[8:15]))
    
    #ax[2].fill_between(np.arange(1,ndays),data_95m[8:15],data_95p[8:15],facecolor=color_stat,alpha=0.15)
    #ax[2].fill_between(np.arange(1,ndays),locm[1:],locp[1:],facecolor=color_stat_spat,alpha=0.15)
    ax[2].set_ylabel('Overlap')
    ax[2].set_xlabel(r'$\Delta t$ (sessions)',labelpad=padx)
    #ax[2].plot(np.arange(1,ndays),avg_chance[1:],'--',linewidth=stat_width,color='gray')
    
    ax[2].set_ylim([0.3,0.9])
    ax[2].set_yticks([0.3,0.9])
    ax[2].set_yticklabels([0.3,0.9])
    ax[2].set_xticks([1,7])
    ax[2].set_xticklabels([1,7])
    if legend:
        ax[1].legend(frameon=False,fontsize=6,loc=[0.25,0.44])
    ##inset #################
    
    axins = ax[2].inset_axes([0.6,0.7,0.4,0.4])
    axins.errorbar(np.arange(0,ndays),data_mice[:,0],yerr=data_mice[:,1],fmt='-o',color=color_data,markersize=2,mfc=scatt_data_col,linewidth=w_line,label=' exp. data')
    
    #local_overl=np.concatenate(([1],avg_solutions_spat[8:15]))
    axins.plot(np.arange(0,ndays),local_overl[:],color=color_stat_spat,label='stat. model',zorder=3,linewidth=stat_width)
    netw_data_overlap=np.row_stack(([1,0],network_data[2*ndays:3*ndays-1,:]))
    axins.errorbar(np.arange(0,ndays),netw_data_overlap[:,0],yerr=netw_data_overlap[:,1],fmt='-*',color=color_network,markersize=2,mfc=scatt_net_col,linewidth=w_line,label='network model')

    #locm=np.concatenate(([1],data_95m_spat[8:15]))
    #locp=np.concatenate(([1],data_95p_spat[8:15]))
    
    #ax[2].fill_between(np.arange(1,ndays),data_95m[8:15],data_95p[8:15],facecolor=color_stat,alpha=0.15)
    #axins.fill_between(np.arange(0,ndays),locm[:],locp[:],facecolor=color_stat_spat,alpha=0.15)
    #axins.plot(np.arange(0,ndays),avg_chance[:],'--',linewidth=1,color='gray')
    
    #axins.set_ylabel('Overlap',fontsize=6)
    #axins.set_xlabel(r'$\Delta t$',labelpad=padx,fontsize=6)
    
    #axins.errorbar(sessions,avg_act_network[:,0],yerr=avg_act_network[:,1],fmt='-*',color=color_network,markersize=2,mfc=scatt_net_col,linewidth=w_line)
    # #axins.plot(sessions,model_avg_act_vec,'-',color=color_stat)
    axins.spines['right'].set_visible(False)
    axins.spines['top'].set_visible(False)
    
    
    axins.set_ylim([0.,1.1])
    axins.set_yticks([0.,1])
    
    axins.set_yticklabels([0.,1],fontsize=6)
    daysnum=[0,7]
    axins.set_xticks(daysnum)
    axins.set_xticklabels(daysnum,fontsize=6)
    axins.set_xlabel(r'$\Delta t$',labelpad=-5,fontsize=6)
    axins.set_ylabel("Overlap ",fontsize=6,labelpad=0)
    axins.tick_params('x',pad=1,length=2)
    axins.tick_params('y',pad=1,length=2)
    
    ###################################################################
    
    
    #ax[2,0].plot(np.arange(1,ndays),avg_over2[1:])
    #ax[2,0].plot(np.arange(1,ndays),avg_over3[1:],color='orange',linewidth=1)
    
    #ax[2,0].set_ylim([0.4,1])
    #ax[2,0].set_yticks([0.4,1])
    
    
    if plotting_exp:    
        ax[0].errorbar(np.arange(1,ndays+1),data_total[2*ndays-1:3*ndays-1,0],yerr=data_total[2*ndays-1:3*ndays-1,1],fmt='-o',color=color_data,label=' exp. data',zorder=1,markersize=ms,mfc=scatt_data_col,linewidth=w_line)
    ax[0].errorbar(np.arange(1,ndays+1),network_data[:ndays,0],yerr=network_data[:ndays,1],fmt='-*',color=color_network,label='network model',zorder=2,markersize=ms2,mfc=scatt_net_col,linewidth=w_line)
    
    #ax[0,0].errorbar(np.arange(1,ndays+1),network_data_sym[:ndays,0],yerr=network_data_sym[:ndays,1],fmt='-*',color=color_network_sym,markersize=ms2,mfc=scatt_net_col_sym,linewidth=w_line_sh,linestyle=st_sh,alpha=alpha_sh)
    
    
    
    #ax[0].fill_between(np.arange(1,ndays+1),data_95m[15:],data_95p[15:],facecolor=color_stat,alpha=0.15)
    #ax[0].fill_between(np.arange(1,ndays+1),data_95m_spat[15:],data_95p_spat[15:],facecolor=color_stat_spat,alpha=0.15)
    
    #ax[0,0].plot(np.arange(1,ndays+1),avg_hist2)
    #ax[0,0].plot(np.arange(1,ndays+1),avg_hist3,color='orange')
    
    ax[0].set_ylabel('Cells (\%)')
    ax[0].set_xlabel('\# Sessions active',labelpad=padx)
    ax[0].set_ylim([0,0.3])
    ax[0].set_yticks([0,0.3])
    ax[0].set_yticklabels([0,30])
    #ax[2].legend(frameon=False,fontsize=6,loc=[0.3,0.35])
    ax[0].set_xticks(xlab)
    
    ax[0].set_xticklabels(xlab)
    
    
    
    #ax[0,1].set_title("Place cells",y=1.3)
    
    
    plt.subplots_adjust(top=0.93,hspace=0.5,right=0.98,wspace=0.4,bottom=0.15,left=0.3)

    #plt.legend()
    #axins.get_legend().remove()
    #plt.tight_layout()
    figure_title1=figure_name+'.png'
    figure_title2=figure_name+'.pdf'
    figure_title3=figure_name+'.svg'

    plt.savefig(figure_title1,format='png',dpi=300)
    plt.savefig(figure_title2,format='pdf',dpi=300)
    
    plt.savefig(figure_title3,format='svg')
    
    plt.close()	