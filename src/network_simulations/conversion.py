#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 15:47:16 2021

@author: federico
"""

import numpy as np
import scipy.io


def raw_data_conv(data,deltat,tend,Nneur):
    tstart=deltat #in milliseconds
    time=np.arange(tstart,tend+deltat,deltat)
    raw_activity=np.zeros((time.shape[0],Nneur))
    for i in range(1,Nneur+1):
        data1=data[data[:,0]==i]
        idt=np.round((data1[:,2]-tstart)/deltat).astype(int)
        raw_activity[idt,i-1]=1
    return raw_activity
    

def raw_data_conv_2(data,deltat,tend,Nneur):
    tstart=deltat #in milliseconds
    time=np.arange(tstart,tend+deltat,deltat)
    raw_activity=np.zeros((time.shape[0],Nneur),dtype=np.int32)
    #idt=np.round((data-tstart)/deltat).astype(int)
    #idt[data==0]=0
    for i in range(Nneur):
        idt_loc=np.round((data[data[:,i]>0,i]-tstart)/deltat).astype(int)
        raw_activity[idt_loc,i]=1
    
    return raw_activity
    
    
def downsampling(array,n,ncols,thresh_calcium):
    new_array_length=int(round(array.shape[0]/n))
    #print(new_array_length)
    new_array=np.ones((new_array_length,ncols))*np.nan
    for i in range(1,new_array_length+1):
        new_array[i-1,:]=np.sum(array[(i-1)*n:(i-1)*n+n,:],axis=0)>thresh_calcium
    return new_array

def converting_data(dataneutr,filename,tend,Nneur):        
    #dataneutr=np.load("../phasesnpy.npy")    
    
    ###create bin array over time
    tstart=10 #milliseconda
    deltat=tstart
    
    raw_activity=raw_data_conv_2(dataneutr,deltat,tend,Nneur)
    np.delete(raw_activity,0,0)
    scipy.io.savemat(filename, dict(raw_activity=raw_activity))

    return 


def converting_data2(dataneutr,filename,tend,Nneur):        
    #dataneutr=np.load("../phasesnpy.npy")    
    
    ###create bin array over time
    tstart=10 #milliseconda
    deltat=tstart
    
    raw_activity=raw_data_conv(dataneutr,deltat,tend,Nneur)
    np.delete(raw_activity,0,0)
    scipy.io.savemat(filename, dict(raw_activity=raw_activity))

    return raw_activity

