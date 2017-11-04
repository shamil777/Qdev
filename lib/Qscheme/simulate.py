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

from collections import OrderedDict

import SI

class Scheme():
    def __init__(self,graph=None,file_path=None):
        self.graph = graph
        self.elements = {}
        self.params = None
        
        if( graph != None ):
            self._init_elements()
        
        self.file_path = file_path
        if( self.graph is None ):
            if( self.file_path != None ):
                self.graph_from_file(file_path)
        
        self.refresh_symbols()
        
        
    def refresh_symbols(self):
        self.Hc_sym_phase = None
        self.Hj_sym_phase = None
        self._construct_Hc_sym_phase()
        self._construct_Hj_sym_phase()
        
        self.Hc_sym_cooperN = None
        self.Hj_sym_cooperN = None
        self._construct_Hc_sym_cooperN()
        self._construct_Hj_sym_cooperN()
                
        self.C_matrix_sym = None
        self._construct_Cap_matrix_sym_from_graph()
        
        self._update_params()
        
    def _update_params(self):
        self.params = OrderedDict()
        for element in self.elements.values():
            for parameter,value in element.params.items():
                self.params[parameter] = value
    
    def assign_subscripts_to_nameGroups(self,group_name_list,subscript_list):
        for element in self.elements.values():
            for name_i,group_name in enumerate(group_name_list):
                if( element.group_name == group_name and subscript_list[name_i] != None ):
                    element.change_subscript( subscript_list[name_i] )
                    break
        
        self.refresh_symbols()

    def assign_values_to_parameters(self,parameters_list,values_list):
        for i,parameter in enumerate(parameters_list):
            if( not isinstance(parameter,sympy.Symbol) ):
                parameters_list[i] = sympy.Symbol(parameter)
                
        for element in self.elements.values():
            for parameter in element.params:
                for param_to_change,value in zip(parameters_list,values_list):
                    if( param_to_change == parameter ):
                        element.params[param_to_change] = value
                        break
                    
        self._update_params()

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
        nodes_N = max(self.graph.nodes())
        charge_ops = self._get_charge_ops(nodes_N, cooper_N)
        
        self.Hc_num_cooperN = qp.tensor(*[qp.qzero(2*cooper_N + 1) for i in range(nodes_N)])
        C_inv = np.array( ((self.C_matrix_sym)**-1).subs(self.params) ).astype(np.float64)
        C_inv *= 2*(SI.e)**2/(SI.h*10**9)
        
        for i in range(nodes_N):
            for j in range(nodes_N):
                self.Hc_num_cooperN += charge_ops[i]*C_inv[i,j]*charge_ops[j]
        
        
    def _construct_Hj_sym_cooperN(self):        
        self.Hj_sym_cooperN = 0
        for element in self.elements.values():
            self.Hj_sym_cooperN += element._Hj_symbol_cooperN()
        
    def _construct_Hj_num_cooperN(self,cooper_N,f_ext):
        nodes_N = max(self.graph.nodes())
        raising_ops,lowering_ops = raising_lowering_ops(nodes_N,cooper_N)
    
        self.Hj_num_cooperN = qp.tensor(*[qp.qzero(2*cooper_N + 1) for i in range(nodes_N)])
        for element in self.elements.values():
            if( isinstance(element,nl.SIS_JJ) ):
                node1 = element.node1 - 1
                node2 = element.node2 - 1
                Ej = element.params[element.Ej]
                
                exp = 1
                if( element.flux_enabled == True ):
                    exp = np.exp(1j*2*np.pi*f_ext)
                
                if( element.node1 != 0 and element.node2 != 0 ):
                    self.Hj_num_cooperN += -Ej/2*(raising_ops[node1]*lowering_ops[node2]*exp \
                                                  + lowering_ops[node1]*raising_ops[node2]/exp)
                else:
                    node = node1 + node2 # equals to nonzero node number
                    print(node1,node2)
                    self.Hj_num_cooperN += -Ej/2*(lowering_ops[node]*exp \
                                                  + raising_ops[node]/exp)
                
    
    def _get_charge_ops( self,nodes_N,cooper_N ):
        charge_ops = []
        for i in range(nodes_N):
            to_tensor_list = []
            for j in range(nodes_N):
                if( j != i ):
                    to_tensor_list.append( qp.identity(2*cooper_N+1) )
                else:
                    to_tensor_list.append( qp.charge(cooper_N) )
            charge_ops.append( qp.tensor(*to_tensor_list) )
        return charge_ops
    
    def _construct_Cap_matrix_sym_from_graph(self):
        ''' 
        @description:
            
        @parameters:
            
        @return:
        '''
        N = len(self.graph)-1
        C_matrix = sympy.zeros( N,N )
        for edge in self.graph.edges(data="element"):
            # to make nodei mutable type
            node1 = edge[0]
            node2 = edge[1] 
            
            element = edge[2]
            if( isinstance(element,(nl.Cap,nl.SIS_JJ)) ):
                C = element.C
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
        self._init_elements()
        

e = 1.6e-19 # electron charge SI
hbar = 1.055e-34 # plank constant SI (this is \hbar)
h = hbar*(2*np.pi)

alpha = 0.43 # ratio of the smaller and larger josephson junctions areas
E_C = 17.6 # GHz, roughly corresponding to C = 4.5 fF
C = 4*e**2/2/(2*np.pi*E_C*hbar*10**9)
E_J = 36.2/alpha # GHz
Csh = 51e-15/C # in units of C (C equals to larger josephson junctions capacitance)

Z_r = 50 # Ohm - resonator line impendance
freq_r = 6.00 # GHz resonator frequency
omerga_r = freq_r*(2*np.pi)
beta = 0.043 # coupling parameter
#rms single-photon voltage in resonator
V_r0 = freq_r*np.sqrt(2*h*Z_r)

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

def raising(cooper_N):
    return qp.Qobj(np.diag(np.ones(2*cooper_N),-1))

def raising_lowering_ops(nodes_N,cooper_N):
    raising_ops = []
    lowering_ops = []
    for i in range(nodes_N):
        to_tensor_list_r = []
        to_tensor_list_l = []
        for j in range(nodes_N):
            if( j != i ):
                to_tensor_list_r.append( qp.identity(2*cooper_N+1) )
                to_tensor_list_l.append( qp.identity(2*cooper_N+1) )
            else:
                to_tensor_list_r.append( raising(cooper_N) )
                to_tensor_list_l.append( lowering(cooper_N) )
        raising_ops.append( qp.tensor(*to_tensor_list_r) )
        lowering_ops.append( qp.tensor(*to_tensor_list_l) )
        
    return raising_ops,lowering_ops
    
def lowering(N):
    return qp.Qobj(np.diag(np.ones(2*N),1))
# tunneling operator is available in QuTip and his result equals to "rising(N) + lowering(N)"

def H_J( E_J, alpha, phase_ext, N ):
    H1 = 0.5*( np.exp(1j*phase_ext)*qp.tensor(lowering(N),raising(N)) \
                         + np.exp(-1j*phase_ext)*qp.tensor(raising(N),lowering(N)) )
    H2 = qp.tensor(qp.identity(2*N+1), 0.5*qp.tunneling(2*N+1,1))
    H3 = qp.tensor(0.5*qp.tunneling(2*N+1,1), qp.identity(2*N+1))
        
    tunneling_ops = (H1,H2,H3)                  
    H_J = qp.tensor(qp.qzero(2*N+1), qp.qzero(2*N+1))
    E_Js = (-1)*np.array([alpha*E_J,E_J,E_J])
    for idx, tunneling_op in enumerate(tunneling_ops):
        H_J += E_Js[idx]*tunneling_op
    return H_J

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