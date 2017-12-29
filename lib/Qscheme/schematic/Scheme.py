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

import itertools


import Qscheme.schematic._netlist_parser as nl
import Qscheme.variables as vs

from collections import OrderedDict

import SI

class Scheme():
    def __init__(self,graph=None,file_path=None):
        self.file_path = file_path
        self.graph = graph
        self.elements = OrderedDict()
        self.params = OrderedDict()
        
        #### SYMBOLS AND VARIABLES DEFINITIONS START ####
        self.Hc_sym_phase = None
        self.Hj_sym_phase = None
        self.Hc_sym_cooperN = None
        self.Hj_sym_cooperN = None
        self.C_matrix_sym = None
        
        self.nodes_marks_list = ["GND"] # first node is always ground node
        self.cooper_N = None
        
        # constructing graph from file, if it is not 
        # passed as an argument
        if( self.graph is None ):
                self.graph_from_file(file_path)
        
        # initializing scheme node variables (only symbolic part)
        self.nodes_N = len(self.graph)        
        self.r_ops, self.l_ops = \
            vs._get_raising_lowering_ops_as_vars(self.nodes_N)
        self.n = vs._get_cooperN_ops_as_vars(self.nodes_N)
        self.phi = vs._get_phase_ops_as_vars(self.nodes_N)
        self.phi_d = vs._phase_pd_ops_as_vars(self.nodes_N)            
        #### SYMBOLS AND VARIABLES DEFINITIONS END ####
        
        # initializing elements dict
        self._init_elements()
        
        
        # Copies all parameters to the self.params dict.
        # Params with the same name will be overwritten.
        # But this is not the case in this moment in constructor
        self._load_params_to_scheme_class()
        
        # supply Elements with reference to this Scheme class
        self._provide_elements_with_scheme_instance()
        
        # clears and constructs new Hc, Hj and Cmatrix symbols
        # based on corresponding graph elements functions
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
        self.H_cooperN = self.Hc_sym_cooperN + self.Hj_sym_cooperN
                
        self.C_matrix_sym = None
        self._construct_Cap_matrix_sym_from_graph()
        
        self.C_matrix_internal_sym = None
    
    
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
    
            
    def _load_params_to_scheme_class(self):
        self.params = OrderedDict()
        for element_sym,element in self.elements.items():
            for elem_param in element.params:
                self.params[elem_param.sym] = elem_param
    
    
    def _init_elements(self):
        '''
        @description:
            initializing scheme.elements dictionary with refDes keys
            and Element class of elements as values
        '''
        for edge in self.graph.edges(data="element"):
            element = edge[2]
            self.elements[element.refDes] = element
    
    
    def _connect_element_params_to_scheme_params(self):
        '''
        @description:
            Storing parameters as classes allows to decrease number
        of parameter classes. Function forces all variables in Element 
        classes to refer to the Scheme.params variables.
            This made for rapid parameter sweep. The parameters 
        one only need to be changed via Scheme.params. 
            So, there is no need to iterate over all graph for each time you 
        need to set new parameter value.
        '''
        for element in self.elements.values():
            for i,param in enumerate(element.params):
                element.params[i] = self.params[param.sym]
            element._connect_to_params_list()

    
    def _provide_elements_with_scheme_instance(self):
        for refDes in self.elements:
            self.elements[refDes].scheme = self
    
    
    def _construct_ops(self,cooper_N):
        '''
        @description:
            constructs raising, lowering and cooperN ops
            provided by cooper_N dimensions
        '''
        if( self.cooper_N == cooper_N  ):
            return
        else:
            self.cooper_N = cooper_N
            self.r_ops,self.l_ops = vs._get_raising_lowering_ops_as_vars(self.nodes_N,self.cooper_N)
            self.n = vs._get_cooperN_ops_as_vars(self.nodes_N,self.cooper_N)
            
    
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
        
        
    def _mark_nodes(self):
        pass
    
    def get_params_values(self):
        return [var.val for var in self.params.values()]
    
    def set_params_values(self, values):
        for i,val in enumerate(values):
            self.params[i] = val