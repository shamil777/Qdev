import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

import Qscheme.simulate.simulate as q

import sympy
import numpy as np

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
    
    cooper_N = 6 
    
    current_dir = "C:\\Users\\user\\Documents\\Qdev\\qubit_simulations\\Qutip_Plus_DipTrace"
    scheme = q.Scheme(file_path=current_dir + "\\" + "3JJ_V.net")
    alpha_J = scheme.elements["J4"]
    J_N = 4
    Jtest = scheme.elements["J" + str(J_N)]
    #battery = scheme.elements["B1"]
    alpha_J.enable_ext_flux()
    
    
    
    group_names = ["JosType1","JosType2","C1"]
    subscripts = ["J\\alpha","J","sh"]
    scheme.assign_subscripts_to_nameGroups(group_names,subscripts) # values refresh
    
    

    parameters = ["C_{J}","C_{J\\alpha}","C_{sh}","E_{J}","E_{J\\alpha}","\\varphi_{J4}"]
    param_vals = [C_Jb,C_Ja,Csh,E_J,alpha*E_J,f_ext]    
    scheme.assign_values_to_parameters(parameters,param_vals)
    
    display(scheme.Hc_sym_cooperN)
    display(scheme.Hj_sym_cooperN)
    display(scheme.C_matrix_sym)
    N_pts = 21


    d1 = []
    f_ext_list = np.linspace(0.4,0.6,N_pts)
    for i,f_ext in enumerate(f_ext_list):
        param_vals[-1] = f_ext
        scheme.assign_values_to_parameters(parameters,param_vals)
        Evs,Evals = scheme.find_eigensystem(cooper_N,4)
        
        d1.append(Evals)
    d1 = np.array(d1)
    q.plot_eigenergies_custom(d1,f_ext_list)
    '''
    '''
    #display(scheme.Hc_sym_cooperN)
    #display(scheme.Hj_sym_cooperN)
    #display(scheme.C_matrix_sym[0])  
    #display(sympy.sympify((scheme.C_matrix_sym)**(-1)))
        