import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

from collections import OrderedDict

from qscheme.simulator import SchemeSimulator
from qscheme.Schematic import Scheme

import numpy as np

from sympy import init_printing
init_printing()

from qscheme.variables import Var

from IPython.display import display

import matplotlib.pyplot as plt

from SI import e,h

if( __name__ == "__main__" ):

    alpha = 0.5
    E_C = 2.5 # GHz
    E_J = 70 # GHz
    f_ext = 0.5
    C_Jb = 4*e.val**2/2/(h.val*E_C*10**9)
    C_Ja = alpha*C_Jb
    Csh = 5e-15
    
    current_dir = "C:\\Users\\botan\\Documents\\Qdev\\qubit_simulations\\Qutip_Plus_DipTrace"
    scheme = Scheme(file_path=current_dir + "\\" + "2_3JJ_Ccoupled_noV.net")
    
    # changing symbols for those elements that belong to the same group
    group_names = ["JJ_mini","JJ_large","C1"]
    subscripts = ["J\\alpha","J","sh"]
    scheme.assign_subscripts_to_nameGroups(group_names,subscripts) # values refresh
    
    # assigning numerical values to the parameters
    parameters = ["C_{J}","C_{J\\alpha}","C_{sh}","E_{J}","E_{J\\alpha}","\Phi_{J3}","\Phi_{J4}","\Phi_{J5}"]
    param_vals = [C_Jb,C_Ja,Csh,E_J,alpha*E_J,f_ext,f_ext,f_ext]
    scheme.assign_values_to_parameters(parameters,param_vals)
    
    ## printing hamiltonian
    #display(scheme.Hc_sym_cooperN)
    #display(scheme.Hj_sym_cooperN)
    #display(scheme.C_matrix_sym)

    ### SIMULATION # ##
    cooper_N = 5
    N_pts = 21
    sim = SchemeSimulator(scheme,cooper_N)
    J3_var = scheme.elements["J3"].f_ext
    J4_var = scheme.elements["J4"].f_ext
    J5_var = scheme.elements["J5"].f_ext
    
    t = Var("t",0)
    scheme.params[t.sym] = t
    mesh = OrderedDict([[J3_var.sym,-np.linspace(0.4,0.6,N_pts)+0.05],
                         [J4_var.sym,0*np.linspace(0.4,0.6,N_pts)],
                         [J5_var.sym,np.linspace(0.4,0.6,N_pts)+0.05]])
    print(scheme.fundamental_cycles_info_str())
    sim.find_eigensystem_internal_parametric(mesh,4)
    
    ### VISUALIZATION ### 
    
    pts = np.array(sim.points)
    #x = pts[:,list(scheme.params.keys()).index(J4_var.sym)]
    x = np.linspace(0.4,0.6,N_pts)
    spectr_idxs = [[0,1],[0,2],[0,3]]

    sim.einvects_list = np.array(sim.einvects_list)
    sim.evals_list = np.array(sim.evals_list)
    y_list = []*len(spectr_idxs)
    for idx in spectr_idxs:
        # y_list contains E_idx - E_0 for all idx in specr_idx
        y_list.append(sim.evals_list[:,idx[1]] - sim.evals_list[:,idx[0]])
    sim.einvects_list = list(sim.einvects_list)
    sim.evals_list = list(sim.evals_list)
    
    # plotting
    for idx,y in zip(spectr_idxs,y_list):
        plt.plot(x,y, label="$E_" + str(idx[1]) + " - E_"+ str(idx[0]) + "$")
        
    # plot cosmetics
    plt.ylabel(r"E, GHz")
    plt.xlabel(r"$" + str(J3_var.sym) + "$")
    plt.legend()
    plt.grid()
    
    plt.show()  

        