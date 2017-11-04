import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

import Qscheme.simulate as q
import numpy as np

import sympy
from sympy import init_printing
init_printing()

from IPython.display import display

from SI import e,hbar,h

if( __name__ == "__main__" ):
    
    C_Jb = 4.5e-15 # F
    alpha = 0.43
    C_Ja = alpha*C_Jb
    Csh = 10*C_Jb
    
    E_C = 4*e**2/(2*C_Jb)/(h*10**9) # GHz
    E_J = 86 # GHz
    
    current_dir = "C:\\Users\\user\\Documents\\Qdev\\qubit_simulations\\Qutip_Plus_DipTrace"
    scheme = q.Scheme(file_path=current_dir + "\\" + "PAD.net")
    
    group_names = ["JosType1","JosType2","C1"]
    subscripts = ["J\\alpha","J",None]
    
    scheme.elements["J4"].enable_ext_flux()
    scheme.assign_subscripts_to_nameGroups(group_names,subscripts) # values refresh
    
    parameters = ["C_{J}","C_{J\\alpha}","C_{12}","E_{J}","E_{J\\alpha}"]
    param_vals = [C_Jb,C_Ja,Csh,E_J,alpha*E_J]
    scheme.assign_values_to_parameters(parameters,param_vals)
    
    scheme._construct_Hc_num_cooperN(1)
    scheme._construct_Hj_num_cooperN(1,0.5)
    print(scheme.Hj_num_cooperN)
    
    #TODO parameters value substitution to produce numberical Hamiltonian

    display(scheme.Hc_sym_cooperN)
    display(scheme.Hj_sym_cooperN)
    display(scheme.C_matrix_sym[0])  
    #display(sympy.sympify((scheme.C_matrix_sym)**(-1)))
        