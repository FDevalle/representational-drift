#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 14:44:23 2021

@author: federico
"""

#import network_sim
import sys
import numpy as np
import time
#import matlab.engine
#import fig_paper
import pathlib


def main_network_simulations(BASE_DIR):
    from src.network_simulations import beta, matrices, conversion
    Nexc = 4000
    Ninh = 1000
    Nff = 8000  # Ca3 cells
    Ncontext = 8000  # EC3 cells
    ca3_spars = 0.5
    frac_tun = 0.5
    frac_tun_ec = 0.
    cont_spars = 0.5  # ec3 sparseness
    nlaps = 1  # number of laps on the track for each session
    ttot = 7.*nlaps*8*1000  # total simulation time
    length_ts = int(100*ttot/1000)
    Nff_act = int(ca3_spars*Nff)
    Ncontext_act = int(cont_spars*Ncontext)
    spatial_neurons = int(frac_tun*Nff_act)
    spatial_neurons_ec = int(frac_tun_ec*Ncontext_act)
    # connectivy matrices
    rhos = 0.95  # ca3 time correlation
    rhons = 0.3  # ec3 time correlation
    Kinh = 500
    fexc = 0.125
    finh = 0.5
    aei, aii, aie, matr_space_inh, matr_cont_inh, nonspatial_matrix, spatial_matrix = matrices.connectivity_matrices(
        Nexc, Nff, Ninh, Ncontext, rhos, rhons, fexc, finh, Kinh)

    np.save(BASE_DIR + '/data/outputs_simulation/aei.npy', aei)
    np.save(BASE_DIR + '/data/outputs_simulation/aie.npy', aie)
    np.save(BASE_DIR + '/data/outputs_simulation/spatial_matrix.npy', spatial_matrix)
    np.save(BASE_DIR + '/data/outputs_simulation/nonspatial_matrix.npy',
            nonspatial_matrix)
    np.save(BASE_DIR + '/data/outputs_simulation/matr_space_inh.npy', matr_space_inh)
    # network simulation
    output_matrix, output_matrix_inh = beta.network_sim(ttot, nlaps, length_ts, Nff_act, spatial_neurons, spatial_neurons_ec, Ncontext_act,
                                                        spatial_matrix, nonspatial_matrix, aei, aii, aie, matr_cont_inh, matr_space_inh, BASE_DIR,
                                                        Nexc, Nff, Ncontext, Ninh)
    np.save(BASE_DIR + '/data/outputs_simulation/output_matrix', output_matrix)
    conversion.converting_data(
        output_matrix, BASE_DIR + '/data/outputs_simulation/time_series.mat', ttot, Nexc)


if __name__ == "__main__":
    import os
    # add root to sys.path
    parents = filter(lambda path: path.name == 'representational-drift',
                     pathlib.Path(__file__).absolute().parents
                     )
    BASE_DIR = next(parents)
    sys.path.append(str(BASE_DIR))
    # compile fortran routine
    net_sim_exit = os.system(
        "f2py3 --debug-capi -c beta.f90 -m beta")
    print("fortan code compiled with exit code ", net_sim_exit)
    print("running simulation")
    start_time = time.time()
    # run simulation
    main_network_simulations(str(BASE_DIR))
    print("simulation time: ", time.time()-start_time, ' seconds')
