#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 15:37:55 2022

@author: federico
"""


import numpy as np 
import matlab.engine 
import matplotlib.pyplot as plt 
import centroid_figure
import plot_pv_correlation
import stat_model_simulations
import plot_statistics_fun
import heat
from scipy.optimize import curve_fit
import set_figure
import plot_ensemble_corr
import plot_tuning_correlation
def double_exponential(t,tau1,a1,tau2,a2):
    return a1*np.exp(-t/tau1)+a2*np.exp(-t/tau2)
def single_exponential(t,tau1,a1):
    return a1*np.exp(-t/tau1)
#modificato per prova
eng = matlab.engine.start_matlab()
folder_names=["/home/federico/postdoc/hippocampus/simulations_paper/simulations_cluster_rhos_values/rhos_0_85/",
              "/home/federico/postdoc/hippocampus/simulations_paper/simulations_cluster_rhos_values/rhos_0_90/",
              "/home/federico/postdoc/hippocampus/simulations_paper/simulations_cluster_rhos_values/rhos_0_95/",
              "/home/federico/postdoc/hippocampus/simulations_paper/simulations_cluster_rhos_values/rhos_0_95_high_inh/",
              "/home/federico/postdoc/hippocampus/simulations_paper/simulations_cluster_rhos_values/rhos_0_95_rhons_0_3/",
              "/home/federico/postdoc/hippocampus/simulations_paper/simulations_cluster_rhos_values/rhos_0_95_rhons_0_35/",
              "/home/federico/postdoc/hippocampus/simulations_paper/simulations_cluster_rhos_values/rhos_0_95_rhons_0_4/",
              "/home/federico/postdoc/hippocampus/simulations_paper/simulations_cluster_rhos_values/rhos_0_95_rhons_0_45/",
              "/home/federico/postdoc/hippocampus/simulations_paper/simulations_cluster_rhos_values/rhos_0_95_rhons_0_5/"]

pv_correlations_file_names=[s + "pv_correlation_within" for s in folder_names]
centroids_file_names=[s + "centroid" for s in folder_names]
heatmap_file_names=[s + "heatmap" for s in folder_names]
time_series_file_names=[s + "time_series.mat" for s in folder_names]
active_statistics_name=[s + "active_stat" for s in folder_names]
ens_correlations_file_names=[s + "ens_correlation_within" for s in folder_names]
tuning_correlations_file_names=[s + "tuning_correlation_within" for s in folder_names]

parameters_stat_model=[[0.507,0.0310,1.77,0.85],
                        [0.894,0.238,1.825,0.9],
                        [0.579,0.52,1.04,0.95],
                        [0.579,0.52,1.04,0.95],
                        [0.73387805,0.3, 1.38486114, 0.95        ],
                        [0.55590571, 0.35,1.16646393, 0.95],
                        [0.51246715,  0.4       ,1.13121065, 0.95      ],
                       [0.59916718, 0.45      ,1.1476384 , 0.95      ] ,
                       [0.6486823 , 0.5       , 1.13386127, 0.95      ]]

optimizing=[False,False,False,False,False,False,False,False,False]
dataexp=np.loadtxt('./data_mice_c16_c6_averaged.dat')
nsess=8
active_mice=dataexp[nsess-1:2*nsess-1]

thresholds=[0.099547,0.10431,0.25156,0.15131,0.17588,0.16837,0.1315,0.12509,0.13236]
pv_corr_matrix=np.zeros((nsess,len(pv_correlations_file_names)))
net_data_act_matrix=np.zeros((31,len(pv_correlations_file_names)))
density_matrix=np.zeros((7,21,len(pv_correlations_file_names)))
cumul_matrix=np.zeros((7,11,len(pv_correlations_file_names)))
density_com_matrix=np.zeros((7,21,len(pv_correlations_file_names)))
pv_tun_corr_matrix=np.zeros((nsess-1,nsess-1,len(pv_correlations_file_names)))
pv_tun_corr_matrix_place=np.zeros((nsess-1,nsess-1,len(pv_correlations_file_names)))
ens_corr_matrix=np.zeros((nsess-1,nsess-1,len(pv_correlations_file_names)))
ens_corr_matrix_place=np.zeros((nsess-1,nsess-1,len(pv_correlations_file_names)))

for i,ts_name in enumerate(time_series_file_names[:]):
    #i=i+4
    [avg_acr_pos_popv,avg_acr_pos_popv_place,density,cumul,net_data_act,
     rate_map_return_1,rate_map_return_8,density_com,
     pv_tun_corr,pv_tun_corr_place,ens_corr,ens_corr_place] = eng.analysis_automatized(
        ts_name,
        optimizing[i],thresholds[i],folder_names[i],False,False,nargout = 12)
    avg_acr_pos_popv=np.squeeze(np.array(avg_acr_pos_popv))
    avg_acr_pos_popv_place=np.squeeze(np.array(avg_acr_pos_popv_place))
    density=np.squeeze(np.array(density))
    cumul=np.squeeze(np.array(cumul))
    net_data_act=np.squeeze(np.array(net_data_act))
    pv_tun_corr=np.squeeze(np.array(pv_tun_corr))
    pv_tun_corr_place=np.squeeze(np.array(pv_tun_corr_place))
    ens_corr=np.squeeze(np.array(ens_corr))
    ens_corr_place=np.squeeze(np.array(ens_corr_place))
    ens_corr_within=np.concatenate(([0.98],np.nanmean(ens_corr,1)))
    ens_corr_within_place=np.concatenate(([0.9744],np.nanmean(ens_corr_place,1)))
    tun_corr_within=np.concatenate(([0.6884],np.nanmean(pv_tun_corr,1)))
    tun_corr_within_place=np.concatenate(([0.8285],np.nanmean(pv_tun_corr_place,1)))

    # centroid_figure.centroid_figure(density,cumul,centroids_file_names[i])
    plot_pv_correlation.plot_pv_correlation(avg_acr_pos_popv, avg_acr_pos_popv_place, 
                                            pv_correlations_file_names[i])
    
    plot_ensemble_corr.plot_ensemble_correlation(ens_corr_within,ens_corr_within_place,
                                                 ens_correlations_file_names[i])
    plot_tuning_correlation.plot_tuning_correlation(tun_corr_within,tun_corr_within_place,
                                                    tuning_correlations_file_names[i])
    stat_model_data=stat_model_simulations.numerical_sim0_plastic_fit(*parameters_stat_model[i])
    model_avg_act_vec_spat=stat_model_data[nsess-1:2*nsess-1]
    net_with_error=np.column_stack((net_data_act,np.zeros(len(net_data_act))))
    plot_statistics_fun.plotting_statistics(dataexp,stat_model_data,net_with_error,
                                            model_avg_act_vec_spat,active_mice,
                                            active_statistics_name[i],ndays=8,plotting_exp=True,
                                        xsize=1.6,ysize=4,legend=False)
    
    pv_corr_matrix[:,i]=avg_acr_pos_popv
    net_data_act_matrix[:,i]=net_data_act
    density_matrix[:,:,i]=density
    cumul_matrix[:,:,i]=cumul
    density_com_matrix[:,:,i]=density_com
    pv_tun_corr_matrix[:,:,i]=pv_tun_corr
    pv_tun_corr_matrix_place[:,:,i]=pv_tun_corr_place
    ens_corr_matrix[:,:,i]=ens_corr
    ens_corr_matrix_place[:,:,i]=ens_corr_place

    data_1=np.reshape(np.array(rate_map_return_1),(20,-1,8),'F')
    data_8=np.reshape(np.array(rate_map_return_8),(20,-1,8),'F')
    # np.save(heatmap_file_names[i]+"_1",data_1)
    # np.save(heatmap_file_names[i]+"_8",data_8)
    # heat.plot_heatmaps(data_1,data_8,folder_names[i]+'heatmap')
# np.save("pv_correlations_matrix",pv_corr_matrix)
# np.save("net_data_act_matrix",net_data_act_matrix)
# np.save("density_matrix",density_matrix)
# np.save("cumul_matrix",cumul_matrix)
# np.save("density_com_matrix",density_com_matrix)
# np.save("ens_corr_matrix",ens_corr_matrix)
# np.save("ens_corr_matrix_place",ens_corr_matrix_place)
# np.save("pv_tun_corr_matrix",pv_tun_corr_matrix)
# np.save("pv_tun_corr_matrix_place",pv_tun_corr_matrix_place)



tuning_correlation=np.nanmean(pv_tun_corr_matrix,axis=1)
tuning_correlation_place=np.nanmean(pv_tun_corr_matrix_place,axis=1)
ens_corr_mean=np.nanmean(ens_corr_matrix,axis=1)
ens_corr_mean_place=np.nanmean(ens_corr_matrix_place,axis=1)

curve_fit(single_exponential,np.arange(1,8),tuning_correlation[:,6])

popt,pcv=curve_fit(double_exponential,np.arange(0,8),ens_corr_matrix[:,6],bounds=([0,100]))


tuncorrfull=np.concatenate(([0.6884],tuning_correlation[:,6]))
tuncorrfullplace=np.concatenate(([0.8284],tuning_correlation_place[:,6]))
enscorrfull=np.concatenate(([0.9884],ens_corr_mean[:,6]))
popttun,pcv=curve_fit(single_exponential,np.arange(0,8),tuncorrfull,bounds=([0,100]))
popttundouble,pcv=curve_fit(double_exponential,np.arange(0,8),tuncorrfull,bounds=([0,100]))

poptcorr,pcv=curve_fit(single_exponential,np.arange(0,8),enscorrfull,bounds=([0,100]))
poptcorrdouble,pcv=curve_fit(double_exponential,np.arange(0,8),enscorrfull,bounds=([0,100]))
#%%
sess=np.arange(8)
single_exp_color="green"
double_exp_color='red'
lw=1
scat_size=3
m_color_data='darkblue'
fig,ax=set_figure.set_figure(2,1,1.5,2.5)
ax[0].plot(sess,enscorrfull,'-o',color=m_color_data,
               linewidth=lw,markersize=scat_size)
##single exponential ##
ax[0].plot(sess,single_exponential(sess,*poptcorr),linewidth=lw,zorder=3,
           color=single_exp_color,label='single exp.')
ax[0].plot(sess,double_exponential(sess,*poptcorrdouble),linewidth=lw,zorder=4,
           color=double_exp_color,label='double exp.')

ax[0].set_xticks([0,7])
ax[0].set_xticklabels([0,7])
ax[0].set_yticks([0.,1])

ax[0].set_yticklabels([0.,1])
ax[0].set_ylabel("Rate correlation")
ax[0].set_xlabel(r"$\Delta$t (sessions)")
ax[0].legend(frameon=False,fontsize=6,loc='upper right')

ax[1].plot(sess,tuncorrfull,'-o',color=m_color_data,
               linewidth=lw,markersize=scat_size)
ax[1].plot(sess,single_exponential(sess,*popttun),linewidth=lw,zorder=3,
           color=single_exp_color,label='single exp.')
ax[1].plot(sess,double_exponential(sess,*popttundouble),linewidth=lw,zorder=4,
           color=double_exp_color,label='double exp.')
ax[1].set_xticks([0,7])
ax[1].set_xticklabels([0,7])
ax[1].set_yticks([0,1])

ax[1].set_yticklabels([0,1])
ax[1].set_ylabel("Tuning correlation")
ax[1].set_xlabel(r"$\Delta$t (sessions)")
ax[1].legend(frameon=False,fontsize=6)

plt.subplots_adjust(left=0.3,bottom=0.15,right=0.95,top=0.95,hspace=1,wspace=0.4)
plt.savefig("curves_exp_network.png",dpi=300)
plt.savefig("curves_exp_network.svg",format='svg')
