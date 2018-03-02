import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

from Qscheme.schematic.Scheme import Scheme
from Qscheme.simulator import SchemeSimulator
from Qscheme.simDialog import SimWindow

from PyQt5 import QtCore,QtWidgets

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
    
    current_dir = "C:\\Users\\botan\\Documents\\Qdev\\qubit_simulations\\Qutip_Plus_DipTrace"
    scheme = Scheme(file_path=current_dir + "\\" + "3JJ_V.net")
    
    cooper_N = 3
    N_pts = 6
    sim = SchemeSimulator(scheme,cooper_N)
    
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()
        
    window = SimWindow(sim)
    sys.exit(app.exec_())

    '''
    ## printing hamiltonian
    display(scheme.Hc_sym_cooperN)
    display(scheme.Hj_sym_cooperN)
    display(scheme.C_matrix_sym)

    ### SIMULATION # ##
    
    print(scheme.fundamental_cycles_info_str())
    print(scheme.variables_info_str())    
    sweep_var = scheme.elements["J4"].f_ext
    mesh = collections.OrderedDict([[sweep_var.sym,np.linspace(0.4,0.6,21)]])
    sim.find_eigensystem_internal_params_product(mesh,4)

    ### VISUALIZATION ###    
    fixed_vars = [var for var in scheme.params.values() if(var.sym != sweep_var.sym)]
    spectr_idxs = [[0,1],[1,2],[2,3]]
    sim.plot2D_evals_from_var(sweep_var,
                              fixed_vars,
                              spectr_idxs)
    
    '''
        