#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 16:55:53 2022

@author: federico
"""
import numpy as np 
import matplotlib.pyplot as plt 
import set_figure

def plot_pv_correlation(pv_sim,pv_sim_place,
                        output_name,ec3_tuned=False,**kwargs):
        
    pv_sim_ec3_tuned = kwargs.get('ec_3_tuned_all_name', None)
    pv_sim_place_tuned = kwargs.get('ec_3_tuned_place_name', None)
    pv_sim[0]=0.86
    pv_sim_place[0]=0.85
    #data_total=np.loadtxt("data_total.dat")
    #active_mice=np.loadtxt("average_active_mice.dat")
    pv_mouse=np.loadtxt("/home/federico/postdoc/Data from Zivlab 2/Data from Zivlab/c16m4/popv_corr_exp_active_cells.dat")
    
    pv_mouse=np.loadtxt('/home/federico/postdoc/Data from Zivlab 2/Data_Liron/drive-download-20210723T135015Z-001/pv_corr_all_cells_AB_all_maps.dat')
    pv_mouse_place=np.loadtxt('/home/federico/postdoc/Data from Zivlab 2/Data_Liron/drive-download-20210723T135015Z-001/pv_corr_place_cells_AB_all_map.dat')
   # pv_mouse=np.loadtxt(filename_data)
    pv_mouse[0,0]=0.67
    pv_mouse_place[0,0]=0.76 ##to be checked
    #pv_mouse[:,0]=1+pv_mouse[:,0]-0.67
    #pv_sim=pv_sim+1-pv_sim[0]
    #pv_mouse=pv_mouse[:,0]
    # pv_sim=np.loadtxt(
    #     "/home/federico/postdoc/hippocampus/plastic_spatial/test_3/new_pars/respecting_distribution/no_shuffling_connections/lower_spat_info/alternative/original_parameters_changing_n_spat/sigma_s_1_91/N_s_1200/not_changing_sigma_d/left_right_running/checking_inhibition/larger_spatial/p_0_13/long_simulation/no_inhibition_heterogeneity/0_5/matlab/PVcorr_simulations_place_cells.dat")
    
    # pv_sim_zero_corr=np.loadtxt("/home/federico/postdoc/hippocampus/plastic_spatial/test_3/new_pars/respecting_distribution/no_shuffling_connections/lower_spat_info/alternative/original_parameters_changing_n_spat/sigma_s_1_91/N_s_1200/not_changing_sigma_d/left_right_running/checking_inhibition/larger_spatial/100_percent/symmetric_solution/matlab/PV_corr_random_sim.dat")
    
    # pvsim_2=np.loadtxt("/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/matlab/PVcorr_simulations_all.dat")
    
    # pv_sim_old=np.loadtxt('/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_rhos_0_95/matlab/pv_corr_all_cells.dat',delimiter=',')
    # pv_sim_place_old=np.loadtxt('/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_rhos_0_95/matlab/pv_corr_place_cells.dat',delimiter=',')
    
    # #pv_sim=np.loadtxt('/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/considering_inhibition/finer_time_resolution/save_inhibitory_currents/matlab/pv_corr_all_cells.dat',delimiter=',')
    # pv_sim_place=np.loadtxt('/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/considering_inhibition/finer_time_resolution/save_inhibitory_currents/matlab/pv_corr_place_cells.dat',delimiter=',')
    # pv_sim_place=np.loadtxt(filename_no_tun_place,delimiter=',')
    # pv_sim_place_tuned=np.loadtxt(filename_tuned_place,delimiter=',')

    
    # pv_sim=np.loadtxt(filename_no_tun,delimiter=',')
    # pv_sim_ec3_tuned=np.loadtxt(filename_tuned,delimiter=',')
    
    
    
    
    
    # ens_corr_sim_all=np.loadtxt('/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_rhos_0_95/matlab/ensemble_corr_all_cells.dat',delimiter=',')
    # ens_corr_sim_place=np.loadtxt('/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_rhos_0_95/matlab/ensemble_corr_place_cells.dat',delimiter=',')
    
    # ens_corr_mouse_all=np.loadtxt('/home/federico/postdoc/Data from Zivlab 2/Data_Liron/drive-download-20210723T135015Z-001/ensemble_corr_all_cells_AB_all_map.dat')
    # ens_corr_mouse_place=np.loadtxt('/home/federico/postdoc/Data from Zivlab 2/Data_Liron/drive-download-20210723T135015Z-001/ensemble_corr_place_cells_AB_all_map.dat')
    
    
    
    # block_CA3_ens=np.loadtxt("/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_simulation_def/blocking_CA3/matlab/ensemble_corr_all_cells.dat",delimiter=',')
    # block_CA3_pv=np.loadtxt("/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_simulation_def/blocking_CA3/matlab/pv_corr_all_cells.dat",delimiter=',')
    
    # block_CA3_ens_place=np.loadtxt("/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_simulation_def/blocking_CA3/matlab/ensemble_corr_place_cells.dat",delimiter=',')
    # block_CA3_pv_place=np.loadtxt("/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_simulation_def/blocking_CA3/matlab/pv_corr_place_cells.dat",delimiter=',')
    
    
    
    # block_EC_ens=np.loadtxt("/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_simulation_def/blocking_EC/matlab/ensemble_corr_all_cells.dat",delimiter=',')
    # block_EC_pv=np.loadtxt("/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_simulation_def/blocking_EC/matlab/pv_corr_all_cells.dat",delimiter=',')
    
    
    # block_EC_ens_place=np.loadtxt("/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_simulation_def/blocking_EC/matlab/ensemble_corr_place_cells.dat",delimiter=',')
    # block_EC_pv_place=np.loadtxt("/home/federico/postdoc/hippocampus/codes_paper/network/fortran_sim/network_simulation_def/blocking_EC/matlab/pv_corr_place_cells.dat",delimiter=',')
    
    
    
    color_data='darkblue'
    color_network='limegreen'
    color_stat='red'
    color_stat_spat='red'
    
    scatt_net_col='limegreen'
    scatt_data_col='darkblue'
    w_line=1
    ms=3
    ms2=4
    padx=-1.5
    
    fig,ax=set_figure.set_figure(2,1,1.5,2.5)
   
    
    
    ndays=8
    sessions=np.arange(1,ndays+1)
    
    #data_time_full=[1.,         0.70779221, 0.57272727, 0.50194805, 0.45194805, 0.4025974,
    # 0.37077922, 0.35064935]
    ax[0].errorbar(np.arange(0,ndays),pv_mouse_place[:,0],yerr=pv_mouse_place[:,1],
                   fmt='-o',color=color_data,markersize=ms,mfc=scatt_data_col,linewidth=w_line,capsize=0.5)
    ax[0].plot(np.arange(0,ndays),pv_sim_place,'-*',color=color_network,markersize=ms2,
               mfc=scatt_net_col,linewidth=w_line,zorder=3,label='control')
    if ec3_tuned:
        ax[0].plot(np.arange(0,ndays),pv_sim_place_tuned,label='tuned ECIII',color='red')
        ax[1].plot(np.arange(0,ndays),pv_sim_ec3_tuned,'-*',color='red',markersize=ms2,mfc='red',linewidth=w_line,zorder=3,label='ecIII tuned')

   # ax[0].plot(np.arange(0,ndays),pv_sim_place_old,'-*',color='orange',markersize=ms2,
               #mfc='orange',linewidth=w_line,zorder=3,label='rho_ns=0.4')
    
    ax[0].set_ylim([0,1.05])
    ax[0].set_yticks([0,1])
    ax[0].set_yticklabels([0,1])
    ax[0].set_ylabel('PV corr')
    ax[0].set_xticks([0,7])
    ax[0].set_xticklabels([])
    ax[0].legend(frameon=False,fontsize=6,loc='lower left')

    ax[1].errorbar(np.arange(0,ndays),pv_mouse[:,0],yerr=pv_mouse[:,1],fmt='-o',color=color_data,markersize=ms,mfc=scatt_data_col,linewidth=w_line,capsize=0.5)
    ax[1].plot(np.arange(0,ndays),pv_sim,'-*',color=color_network,markersize=ms2,mfc=scatt_net_col,linewidth=w_line,zorder=3,label='network')
    #ax[1].plot(np.arange(0,ndays),pv_sim_ec3_tuned,'-*',color='red',markersize=ms2,mfc='red',linewidth=w_line,zorder=3,label='ecIII tuned')
    #ax[1].plot(np.arange(0,ndays),pv_sim_old,'-*',color='orange',markersize=ms2,mfc='orange',linewidth=w_line,zorder=3,label='rho_s=0.4')

    #ax[1].legend(frameon=False)
    
    
    ax[1].set_ylabel('PV corr.')
    ax[1].set_xlabel(r'$\Delta t$ (sessions)',labelpad=padx)
    
    ax[1].set_ylim([np.min(pv_mouse[:,0])-0.1,1.1])
    ax[1].set_yticks([np.round(np.min(pv_mouse[:,0]),decimals=1),1])
    ax[1].set_yticklabels([np.round(np.min(pv_mouse[:,0]),decimals=1),1])
    ax[1].set_xticks([0,7])
    ax[1].set_xticklabels([0,7])
    
    #fig.delaxes(ax[0])
    
    #plt.legend()
    #axins.get_legend().remove()
    #plt.tight_layout()
    plt.savefig(output_name+'.pdf',format='pdf',bbox_inches="tight")
    plt.savefig(output_name+'.svg',format='svg',bbox_inches="tight")

    #plt.savefig("fit_stat_err_data_stat_paper_PV.svg",format='svg')
    
    plt.close()	
    