#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 11:53:35 2020

@author: federico
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

def adding_conn_vm(connections,n_to_add,kappa,Nspat):
    ones=connections==1
    idx=np.arange(0,Nspat)
    idx_ones=idx[ones]
    to_add=np.random.vonmises(0,kappa,n_to_add)
    id_num=np.rint((to_add+np.pi)*Nspat/(2*np.pi)).astype(int)
    to_flip=np.unique(id_num)
    if to_flip.shape[0]<n_to_add:
        cond=(np.unique(to_flip).shape[0])!=n_to_add
        while cond:
            extra=np.random.vonmises(0,kappa,size=n_to_add-to_flip.shape[0])
            extra=np.rint((extra+np.pi)*Nspat/(2*np.pi)).astype(int)
            to_flip=np.concatenate((to_flip,extra))
            to_flip=np.unique(to_flip)
            cond=(np.unique(to_flip).shape[0])!=n_to_add

    inters,x_i,y_i=np.intersect1d(to_flip,idx_ones,return_indices=True)
    while inters.shape[0]>0:
        to_flip=np.delete(to_flip,x_i)
        cond2=True
        ccc=0
        extra=np.zeros(0)
        while cond2:
            to_flip=np.delete(to_flip,np.arange(-1,-extra.shape[0]-1,-1))
            ccc=ccc+1
            #if ccc>1 :
            #print("while",ccc,np.unique(to_flip).shape[0])
            extra=np.random.vonmises(0,kappa,size=inters.shape[0])
            extra=np.rint((extra+np.pi)*Nspat/(2*np.pi)).astype(int)
            to_flip=np.concatenate((to_flip,extra))
            cond2=np.unique(to_flip).shape[0]!=n_to_add


        inters,x_i,y_i=np.intersect1d(to_flip,idx_ones,return_indices=True)

    zeros=connections==0
    idx_zeros=idx[zeros]
    zz=np.setdiff1d(idx_zeros,to_flip)
    to_flip[to_flip==Nspat]=zz[-1]
    #if ((n_to_add-np.unique(to_flip).shape[0])!=0):
    #    print('problem')
    new_connections=np.copy(connections)
    new_connections[to_flip]=1

    return new_connections

def connectivity_matrices(Nca1,Nspat,Nca1inh,Ncontext,rhos,rhons,fexc,finh,Kinh):

    #test_sim=1
    #Nca1=4000
    #Nspat=7692
    #Ncontext=7547
    #Nca1inh=1000
    matr_cont_inh=np.zeros((Nca1inh,Ncontext),dtype=np.int32,order='F')
    rnd_seq=np.random.rand(Nca1inh,Ncontext)
    matr_cont_inh[rnd_seq<fexc]=1
    #np.savetxt("matr_cont_inh.dat",matr_cont_inh,fmt='%i')

    aei=np.zeros((Nca1,Nca1inh),dtype=np.int32,order='F')
    idxinh=np.arange(0,Nca1inh)
    for i in range(1,Nca1+1):

        aei[i-1,np.random.choice(idxinh,size=Kinh,replace=False)]=1
    #rnd_seq=np.random.rand(Nca1,Nca1inh)
    #aei[rnd_seq<500/1000]=1
    #np.savetxt('aei.dat',aei,fmt='%i')


    aie=np.zeros((Nca1inh,Nca1),dtype=np.int32,order='F')
    rnd_seq=np.random.rand(Nca1inh,Nca1)
    aie[rnd_seq<fexc]=1
    #np.savetxt('aie.dat',aie,fmt='%i')


    aii=np.zeros((Nca1inh,Nca1inh),dtype=np.int32,order='F')
    rnd_seq=np.random.rand(Nca1inh,Nca1inh)
    aii[rnd_seq<finh]=1
    #np.savetxt('aii.dat',aii,fmt='%i')

    matr_space_inh=np.zeros((Nca1inh,Nspat),dtype=np.int32,order='F')
    rnd_seq=np.random.rand(Nca1inh,Nspat)
    matr_space_inh[rnd_seq<fexc]=1
    #np.savetxt("ceffinh.dat",matr_space_inh,fmt='%i')



    matr_cont_0=np.zeros((Nca1,Ncontext),dtype=np.int32,order='F')
    p=1.-rhons**2
    alpha=fexc
    rnd_vec=np.random.rand(Nca1,Ncontext)
    matr_cont_0[rnd_vec<alpha]=1
    matr_cont_full=np.zeros((Nca1*8,Ncontext),dtype=np.int32,order='F')
    for i in range(1,9):


        matr_cont_1=np.zeros((Nca1,Ncontext),dtype=np.int32,order='F')

        p11=(1-alpha)*np.sqrt(1-p)+alpha
        p10=alpha*(1-np.sqrt(1-p))

        matr_cont_0flat=matr_cont_0.flatten()
        matr_cont_1flat=matr_cont_1.flatten()

        rnd_seq1=np.random.rand(matr_cont_0flat.shape[0])
        rnd_seq0=np.random.rand(matr_cont_0flat.shape[0])

        matr_cont_1flat[(matr_cont_0flat==1) & (rnd_seq1<p11)]=1
        matr_cont_1flat[(matr_cont_0flat==0) & (rnd_seq0<p10)]=1
        matr_cont_0=matr_cont_1flat.reshape((Nca1,Ncontext))
        #np.savetxt("nonspatial_connectivity%i_%i.dat" %(test_sim-1,i),matr_cont_0,fmt='%i')
        matr_cont_full[(i-1)*Nca1:Nca1*i,:]=matr_cont_0


    matr_cont_0=np.zeros((Nca1,Nspat),dtype=np.int32,order='F')
    p=1.-rhos**2
    alpha=fexc
    rnd_vec=np.random.rand(Nca1,Nspat)
    matr_cont_0[rnd_vec<alpha]=1
    matr_complete=np.zeros((Nca1*Nspat,30))
    matr_spatial_full=np.zeros((Nca1*8,Nspat),dtype=np.int32,order='F')

    for i in range(1,9):

        matr_cont_1=np.zeros((Nca1,Nspat),dtype=np.int32,order='F')

        p11=(1-alpha)*np.sqrt(1-p)+alpha
        p10=alpha*(1-np.sqrt(1-p))

        matr_cont_0flat=matr_cont_0.flatten()
        matr_cont_1flat=matr_cont_1.flatten()

        rnd_seq1=np.random.rand(matr_cont_0flat.shape[0])
        rnd_seq0=np.random.rand(matr_cont_0flat.shape[0])

        matr_cont_1flat[(matr_cont_0flat==1) & (rnd_seq1<p11)]=1
        #matr_cont_1flat[(matr_cont_0flat==1) & (rnd_seq1>=p11)]=0
        matr_cont_1flat[(matr_cont_0flat==0) & (rnd_seq0<p10)]=1

        matr_cont_0=matr_cont_1flat.reshape((Nca1,Nspat))
        matr_complete[:,i-1]=matr_cont_1flat
        matr_spatial_full[(i-1)*Nca1:Nca1*i,:]=matr_cont_0

        #np.savetxt("spatial_connectivity%i_%i.dat" %(test_sim-1,i),matr_cont_0,fmt='%i')

    return aei,aii,aie,matr_space_inh,matr_cont_inh,matr_cont_full,matr_spatial_full