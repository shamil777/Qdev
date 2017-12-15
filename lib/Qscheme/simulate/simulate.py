# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 20:17:14 2017

@author: user
"""
import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

import numpy as np
import sympy
import qutip as qp
import itertools

import matplotlib.pyplot as plt

import Qscheme.netlist as nl
import Qscheme.simulate.operator_constructors as ops
import Qscheme.variables as vs

from collections import OrderedDict

import SI

class Scheme():
    def __init__(self,graph=None,file_path=None):
        self.graph = graph
        self.elements = {}
        
        self.file_path = file_path
        if( self.graph is None ):
                self.graph_from_file(file_path)
        
        self._init_elements()
        
        self.nodes_N = len(self.graph)
        self.cooper_N = None
        self.r_ops, self.l_ops = \
            vs._get_raising_lowering_ops_as_vars(self.nodes_N)
        self.n = vs._get_cooperN_ops_as_vars(self.nodes_N)
        self.phi = vs._get_phase_ops_as_vars(self.nodes_N)
        self.phi_d = vs._phase_pd_ops_as_vars(self.nodes_N)
        
        self.nodes_marks_list = ["GND"] # first node is always ground node
        
        # VARS of the raising, loweriing, node cooper pairs number operator
        # phase and phase derivative operators        
        self._change_symbols_based_on_nodes_marks()
        
        self.params = None
        # Copies all parameters to the self.params dict.
        # Params with the same name will be overwritten.
        # But this is not the case in the moment
        # right after construction.
        self._load_params_to_scheme_class()
        
        # clears and constructs new Hc, Hj and Cmatrix symbols
        # based on corresponding graph elements functions
        self._provide_elements_with_scheme_instance()
        self.refresh_symbols()
        
        
        
    def _mark_nodes(self):
        pass
    
    def _change_symbols_based_on_nodes_marks(self):
        for element in self.elements.values():
            if( isinstance(element,nl.Battery) ):
                n2 = element.node2
                #self.n[n2].sym = sympy.Symbol("V_{" + str(element.node2) + "}")
    
    def _provide_elements_with_scheme_instance(self):
        for refDes in self.elements:
            self.elements[refDes].scheme = self
    
    def _construct_ops(self,cooper_N):
        if( self.cooper_N == cooper_N  ):
            return
        else:
            self.cooper_N = cooper_N
            self.r_ops,self.l_ops = vs._get_raising_lowering_ops_as_vars(self.nodes_N,self.cooper_N)
            self.n = vs._get_cooperN_ops_as_vars(self.nodes_N,self.cooper_N)
    
    def refresh_symbols(self):
        self.Hc_sym_phase = None
        self.Hj_sym_phase = None
        self._construct_Hc_sym_phase()
        self._construct_Hj_sym_phase()
        
        self.Hc_sym_cooperN = None
        self.Hj_sym_cooperN = None
        self._construct_Hc_sym_cooperN()
        self._construct_Hj_sym_cooperN()
        self.H_cooperN = self.Hc_sym_cooperN + self.Hj_sym_cooperN
                
        self.C_matrix_sym = None
        self._construct_Cap_matrix_sym_from_graph()
        
        self.C_matrix_internal_sym = None
        
        
    def _load_params_to_scheme_class(self):
        self.params = OrderedDict()
        for element_sym,element in self.elements.items():
            for elem_param in element.params:
                self.params[elem_param.sym] = elem_param
    
    def _connect_element_params_to_scheme_params(self):
        for element in self.elements.values():
            for i,param in enumerate(element.params):
                element.params[i] = self.params[param.sym]
            element._connect_to_params_list()
    
    def assign_subscripts_to_nameGroups(self,group_name_list,subscript_list):
        for element in self.elements.values():
            for name_i,group_name in enumerate(group_name_list):
                if( element.group_name == group_name and subscript_list[name_i] != None ):
                    element.change_subscript( subscript_list[name_i] )
                    break    
                
        self._load_params_to_scheme_class()
        self._connect_element_params_to_scheme_params()
        self.refresh_symbols()

    def assign_values_to_parameters(self,parameters_list,values_list):
        for i,parameter in enumerate(parameters_list):
            if( not isinstance(parameter,sympy.Symbol) ):
                parameters_list[i] = sympy.Symbol(parameter)
                
        for i,sym_to_assign in enumerate(parameters_list):
            self.params[sym_to_assign].val = values_list[i]

    def _init_elements(self):
        for edge in self.graph.edges(data="element"):
            element = edge[2]
            self.elements[element.refDes] = element
    
    def _construct_Hc_sym_phase(self):
        self.Hc_sym_phase = 0        
        for element in self.elements.values():
            self.Hc_sym_phase += element._Hc_symbol_phase()
        
            
    def _construct_Hc_num_phase(self):
        return NotImplementedError
        
    def _construct_Hj_sym_phase(self):
        self.Hj_sym_phase = 0
        for element in self.elements.values(): 
            self.Hj_sym_phase += element._Hj_symbol_phase()
            
    def _construct_Hj_num_phase(self):
        return NotImplementedError
            
    def _construct_Hc_sym_cooperN(self):
        self.Hc_sym_cooperN = 0
        for element in self.elements.values():
            self.Hc_sym_cooperN += element._Hc_symbol_cooperN()
            
    def _construct_Hc_num_cooperN(self, cooper_N):
        self._construct_ops(cooper_N)
        self.Hc_num_cooperN = 0
        subs_dict = {var.sym:var.val for var in self.params.values()}
        C_inv_num = np.array( ((self.C_matrix_sym)**-1).subs(subs_dict) ).astype(np.float64)
        C_inv_num *= 2*(SI.e.val)**2/(SI.h.val*10**9)   
        for i in range(1,self.nodes_N):
            for j in range(1,self.nodes_N):
                ni = self.n[i].val
                nj = self.n[j].val
                self.Hc_num_cooperN += ni*C_inv_num[i-1,j-1]*nj
          
    def _construct_Hj_sym_cooperN(self):        
        self.Hj_sym_cooperN = 0
        for element in self.elements.values():
            self.Hj_sym_cooperN += element._Hj_symbol_cooperN()
        
    def _construct_Hj_num_cooperN(self,cooper_N):
        self._construct_ops(cooper_N)      
        self.Hj_num_cooperN = 0
        for element in self.elements.values():
            self.Hj_num_cooperN += element._Hj_num_cooperN()                

    
    def _construct_Cap_matrix_sym_from_graph(self):
        ''' 
        @description:
            
        @parameters:
            
        @return:
        '''
        N = self.nodes_N-1
        C_matrix = sympy.zeros( N,N )
        for edge in self.graph.edges(data="element"):
            # to make nodei mutable type
            node1 = edge[0]
            node2 = edge[1] 
            
            element = edge[2]
            if( isinstance(element,(nl.Cap,nl.SIS_JJ)) ):
                C = element.C.sym
                if( node1 != 0 and node2 != 0 ):
                    C_matrix[node1-1,node2-1] -= C
                    C_matrix[node2-1,node1-1] -= C
                if( node1 != 0 ):
                    C_matrix[node1-1,node1-1] += C
                if( node2 != 0 ):
                    C_matrix[node2-1,node2-1] += C

        self.C_matrix_sym = C_matrix
            
        
    
    def graph_from_file(self, file_path):
        self.graph = nl.build_graph_from_file_PADS(file_path)
        self.file_path = file_path
        
    def find_eigensystem(self,cooperN,eigvals_N):
        self._construct_Hc_num_cooperN(cooperN)
        self._construct_Hj_num_cooperN(cooperN)
        H =  self.Hc_num_cooperN + self.Hj_num_cooperN
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
            
        return einvects, np.array(evals)
        
        



def get_inv_dl_caps_matrix( Csh, alpha ):
    return 1/(2*(alpha + Csh) + 1)*np.array( [[1+alpha+Csh,alpha+Csh],[alpha + Csh,1 + alpha + Csh]] )

def H_c( Ec, Csh, alpha, N ):
    islands_N = 2
    charge_ops = qp.tensor(qp.charge(N),
                           qp.identity(2*N+1)), \
                 qp.tensor(qp.identity(2*N+1),
                           qp.charge(N))
    
    H_c = qp.tensor(qp.qzero(2*N+1), qp.qzero(2*N+1))
    inv_dl_caps_matrix = get_inv_dl_caps_matrix( Csh, alpha )
    for i in range(0,islands_N):
        for j in range(0,islands_N):
            H_c += Ec*charge_ops[i]*inv_dl_caps_matrix[i,j]*charge_ops[j]
    return H_c




def array_eins(f_0,f_1,Nf,E_C,E_J,Csh,alpha,cooper_N,eigvals_N):
    arr_evals = []
    arr_einvects = []
    fs = np.linspace(f_0,f_1,Nf)
    for f in fs: 
        print('\r{:g} '.format(f), end='')
        H = H_J(E_J,alpha,f,cooper_N) + H_c(E_C,Csh,alpha,cooper_N)
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
        arr_evals.append(evals)
        arr_einvects.append(einvects)
    return np.array(arr_einvects), np.array(arr_evals), np.array(fs)

def plot_eigenergies(engs,fs):
    for idx in range(1,engs.shape[1]):
        plt.plot(fs,engs[:,idx]-engs[:,0], label=r'$E_{%s}-E_{%s}$'%(idx,0))
    plt.plot( fs,engs[:,2]-engs[:,1], label=r"$E_2 - E_1$" )
    plt.plot( fs,(engs[:,3]-engs[:,0])/2, label=r"$(E_3 - E_0)/2$" )
    plt.plot( fs,(engs[:,2]-engs[:,0])/2, label=r"$(E_2 - E_0)/2$" )
    plt.legend()
    plt.grid()
    
def plot_eigenergies_custom(engs,fs):
    for idx in range(1,3):
        plt.plot(fs,engs[:,idx]-engs[:,0], label=r'$E_{%s}-E_{%s}$'%(idx,0))
    plt.plot( fs,engs[:,2]-engs[:,1], label=r"$E_2 - E_1$" )
    plt.plot( fs,(engs[:,3]-engs[:,0])/2, label=r"$(E_3 - E_0)/2$" )
    plt.plot( fs,(engs[:,2]-engs[:,0])/2, label=r"$(E_2 - E_0)/2$" )
    plt.legend()
    plt.grid()


def array_eins_from_params_product(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,eigvals_N_pts):
    arr_evals = []
    arr_einvects = []
    
    fparams = [fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,eigvals_N_pts]
    for i,param in enumerate(fparams):
        if( not isinstance(param,list) and not isinstance(param,np.ndarray)):
            fparams[i] = [param]
            
    # sorted order is preserved
    prod = itertools.product(*fparams)
    
    for params in prod:
        f = params[0]*2*np.pi
        E_C = params[1]
        E_J = params[2]
        Csh = params[3]
        alpha = params[4]
        cooper_N = params[5]
        eigvals_N = params[6]
        
        H = H_J(E_J,alpha,f,cooper_N) + H_c(E_C,Csh,alpha,cooper_N)
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
        arr_evals.append(evals)
        arr_einvects.append(einvects)
        
    return arr_einvects, np.array(arr_evals)

def array_eins_from_params_lists(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,eigvals_N_pts):
    arr_evals = []
    arr_einvects = []
    # sorted order is preserved
    
    params = [fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,eigvals_N_pts]
    lengths = []
    for i, param in enumerate(params):
        if( isinstance(param,list) ):
            lengths.append(len(param))
        elif( isinstance(param,np.ndarray) ):
            lengths.append(param.shape[0])
        else:
            lengths.append(1)
            
    max_len = max(lengths)
    
    for i,param in enumerate(params):
        if( not isinstance(param,(list, np.ndarray)) ):
            params[i] = itertools.repeat(param,max_len)
        
    iter_list = zip( *params )
    
    for i,vals in enumerate(iter_list):
        f = vals[0]*2*np.pi
        E_C = vals[1]
        E_J = vals[2]
        Csh = vals[3]
        alpha = vals[4]
        cooper_N = vals[5]
        eigvals_N = vals[6]
        
        print('\r{:g} {:g} {:g} {:g} {:g} {:g} {:g}'.format(f,E_C,E_J,Csh,alpha,cooper_N,eigvals_N), end='')
        H = H_J(E_J,alpha,f,cooper_N) + H_c(E_C,Csh,alpha,cooper_N)
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
        arr_evals.append(evals)
        arr_einvects.append(einvects)
        
    return arr_einvects, np.array(arr_evals)

def array_coupling_from_parameters_product(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,multiplier):
    multiplier = 2*e*V_r0*beta/h # GHz
    
    g = []
    
    fparams = [fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,2]
    for i,param in enumerate(fparams):
        if( not isinstance(param,list) and not isinstance(param,np.ndarray)):
            fparams[i] = [param]
            
    # sorted order is preserved
    prod = itertools.product(*fparams)
    
    Evs, Engs = array_eins_from_params_product(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,2)
    Evs = np.array(Evs)
    
    for i,params in enumerate(prod):
        N = params[5]
        n1_op = qp.tensor(qp.charge(N),qp.identity(2*N+1))
        n2_op = qp.tensor(qp.identity(2*N+1),qp.charge(N))
        
        vec0 = Evs[i,0]
        vec1 = Evs[i,1]
        g.append( multiplier*(n1_op-n2_op).matrix_element(vec0.dag(),vec1) )
                 
    return np.array(g)

def array_coupling_from_parameters_lists(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,multiplier):
    multiplier = 2*e*V_r0*beta/h # GHz
    
    g = []
    
    Evs, Engs = array_eins_from_params_lists(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,2)
    Evs = np.array(Evs)
    
    params = [fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,2]
    lengths = []
    for i, param in enumerate(params):
        if( isinstance(param,list) ):
            lengths.append(len(param))
        elif( isinstance(param,np.ndarray) ):
            lengths.append(param.shape[0])
        else:
            lengths.append(1)
            
    max_len = max(lengths)
    
    for i,param in enumerate(params):
        if( not isinstance(param,(list, np.ndarray)) ):
            params[i] = itertools.repeat(param,max_len)
        
    iter_list = zip( *params )
    
    for i, params in enumerate(iter_list):
        N = params[5]
        n1_op = qp.tensor(qp.charge(N),qp.identity(2*N+1))
        n2_op = qp.tensor(qp.identity(2*N+1),qp.charge(N))
        
        vec0 = Evs[i,0]
        vec1 = Evs[i,1]
        g.append( multiplier*(n1_op-n2_op).matrix_element(vec0.dag(),vec1) )
                 
    return np.array(g)

def array_eins_full_from_parameters_product(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts):
    pass