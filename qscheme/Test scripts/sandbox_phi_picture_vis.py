# -*- coding: utf-8 -*-
import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

import qscheme.simulator.SchemeSimulator as qss
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from math import atan2



if( __name__ == "__main__" ):
    alpha = 0.8 # ratio of the smaller and larger josephson junctions areas
    E_C = 17.6  # GHz, roughly corresponding to C = 4.5 fF
    E_J = 86 # GHz
    Csh = 10.0 # in units of C (C equals to larger josephson junctions capacitance)
    f = 0.5 # flux
    cooper_N = 15
    N_funcs = 2
    Evs, Engs = qss.array_eins_from_params_product(f,E_C,E_J,Csh,alpha,cooper_N,N_funcs)
    N_f = 100
    psi_list = [np.zeros( (N_f,N_f),dtype=np.complex128 ) for i in range(N_funcs)]
    
    f_list = np.linspace( -np.pi,np.pi, N_f )
    f1,f2 = np.meshgrid(f_list,f_list)
    
    evs_list = [np.reshape(Evs[0][i].data.todense(),(2*cooper_N+1,2*cooper_N+1)) for i in range(0,N_funcs)]

    for i in range(N_funcs):
        evs = evs_list[i]
        for n1 in range(-cooper_N,cooper_N+1):
            for n2 in range(-cooper_N,cooper_N+1):
                psi_list[i] += evs[n1+cooper_N,n2+cooper_N]*np.exp(1j*n1*f1 + 1j*n2*f2)
    
    fig = plt.figure()
    rows_max = 3
    cols_max = 3
    cols = cols_max
    if( N_funcs < cols_max ):
        cols = N_funcs    
    
    axes = []
    lines = []
    for i in range(N_funcs):
        axes.append(fig.add_subplot(int(N_funcs/cols_max)+1,cols,i+1,projection="3d"))
        lines.append(axes[i].plot_surface(f1,f2,np.abs(psi_list[i])))
    fig.show()
    