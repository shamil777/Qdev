import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

from Qscheme.schematic.Scheme import Scheme
import Qscheme.simulator.SchemeSimulator as ss

import matplotlib.pyplot as plt
import numpy as np

import collections

from sympy import init_printing
init_printing()

from IPython.display import display


from SI import e,h

if( __name__ == "__main__" ):

    alpha = 0.43 
    E_C = 17.6 # GHz
    E_J = 86 # GHz
    f_ext = 0.5    
    C_Jb = 4*e.val**2/2/(h.val*E_C*10**9)
    C_Ja = alpha*C_Jb
    Csh = 51e-15
    
    current_dir = "C:\\Users\\user\\Documents\\Qdev\\qubit_simulations\\Qutip_Plus_DipTrace"
    scheme = Scheme(file_path=current_dir + "\\" + "3JJ_noV.net")
    alpha_J = scheme.elements["J4"]
    
    # changing symbols for those elements that belong to the same group
    group_names = ["JosType1","JosType2","C1"]
    subscripts = ["J\\alpha","J","sh"]
    scheme.assign_subscripts_to_nameGroups(group_names,subscripts) # values refresh
    
    # assigning numerical values to the parameters
    parameters = ["C_{J}","C_{J\\alpha}","C_{sh}","E_{J}","E_{J\\alpha}","\Phi_{J4}"]
    param_vals = [C_Jb,C_Ja,Csh,E_J,alpha*E_J,f_ext]
    scheme.assign_values_to_parameters(parameters,param_vals)
    alpha_J.enable_ext_flux()
    
    # printing hamiltonian
    display(scheme.Hc_sym_cooperN)
    display(scheme.Hj_sym_cooperN)
    display(scheme.C_matrix_sym)

    ### SIMULATION ###
    cooper_N = 4
    N_pts = 21
    sim = ss.SchemeSimulator(scheme,cooper_N)
    
    
    mesh = collections.OrderedDict([["\Phi_{J4}",np.linspace(0.4,0.6,21)]])
    
    sim.find_eigensystem_params_product(mesh,3)
    x = mesh["\Phi_{J4}"]
    evals_list = np.array(sim.evals_list)
    y = evals_list[:,1] - evals_list[:,0]
    plt.plot( x, y )
    plt.show()
    

    #display(scheme.Hc_sym_cooperN)
    #display(scheme.Hj_sym_cooperN)
    #display(scheme.C_matrix_sym[0])  
    #display(sympy.sympify((scheme.C_matrix_sym)**(-1)))
        