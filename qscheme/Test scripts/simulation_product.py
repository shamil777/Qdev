import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

from collections import OrderedDict

from qscheme.simulator import SchemeSimulator
from qscheme.Schematic import Scheme

import numpy as np

from sympy import init_printing
init_printing()

from IPython.display import display

from SI import e,h

if( __name__ == "__main__" ):

    alpha = 0.43 
    E_C = 5 # GHz
    E_J = 70 # GHz
    f_ext = 0.5
    C_Jb = 4*e.val**2/2/(h.val*E_C*10**9)
    C_Ja = alpha*C_Jb
    Csh = 51e-15
    
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
    display(scheme.Hc_sym_cooperN)
    display(scheme.Hj_sym_cooperN)
    display(scheme.C_matrix_sym)

    ### SIMULATION # ##
    cooper_N = 3
    N_pts = 6
    sim = SchemeSimulator(scheme,cooper_N)  
    sweep_var = scheme.elements["J4"].f_ext
    mesh = OrderedDict([[sweep_var.sym,np.linspace(0.4,0.6,21)]])
    sim.find_eigensystem_internal_product(mesh,4)

    ### VISUALIZATION ###    
    fixed_vars = [var for var in scheme.params.values() if(var.sym != sweep_var.sym)]
    spectr_idxs = [[0,1],[1,2],[2,3]]
    sim.plot2D_evals_from_var(sweep_var,
                              fixed_vars,
                              spectr_idxs)

        