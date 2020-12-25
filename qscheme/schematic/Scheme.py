# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 20:17:14 2017

@author: user
"""

from collections import OrderedDict

import numpy as np
import sympy
import networkx as nx

import SI
from qscheme.Schematic import _netlist_parser as nl
import qscheme.variables as vs


class Scheme():
    def __init__(self, graph=None, file_path=None):
        self.file_path = file_path
        self.graph = graph
        self.internal_scheme = None
        
        self.min_span_tree = None
        self.min_span_tree_elements = OrderedDict()
        self.flux_elements = OrderedDict()
        self.fundamental_cycles_elements = OrderedDict()
        self.elements = OrderedDict()
        self.ext_V_elements = OrderedDict()
        self.params = OrderedDict()
        
        # SYMBOLS AND VARIABLES DEFINITIONS START #
        self.Hc_sym_phase = None
        self.Hj_sym_phase = None
        self.Hc_sym_cooperN = None
        self.Hj_sym_cooperN = None
        self.H_coupling_sym_phase = None
        self.H_coupling_sym_cooperN = None
        self.H_external_sym_phase = None
        self.H_external_sym_cooperN = None
        self.C_matrix_sym = None
        self.C_intermediate_sym = None
        self.C_matrix_internal_sym = None
        
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
        
        # finding some spanning tree
        self.min_span_tree = nx.minimum_spanning_tree(self.graph)
        # filling spanning tree elements reference dict
        for edge in self.min_span_tree.edges(data="element"):
            element = edge[2]
            self.min_span_tree_elements[element.refDes] = element
        # filling flux_elements reference dict
        for element_key in self.elements:
            element = self.elements[element_key]
            if( element_key not in self.min_span_tree_elements ):
                self.flux_elements[element.refDes] = element
                
        # finding fundamental cycles that correspond to the spanning tree
        self.find_fundamental_cycles()
        
        # enabling fluxes through all the flux_elements members
        for element_key in self.flux_elements:
            element = self.flux_elements[element_key]
            element.flux_enable()
        
        # Copies all parameters to the self.params dict.
        # Params with the same name will be overwritten.
        # But this is not the case in this moment in constructor
        self._load_params_to_scheme_class()
        
        # clears and constructs new Hc, Hj and Cmatrix symbols
        # also constructing H_coupling and H_external
        # based on corresponding graph elements functions        
        self.refresh_symbols()
        
        # temporary invalid code:
        '''
        temp_graph = None
        for element_class in self.scheme.elements.values():
            if( isinstance(element_class,ExternalForceElement) ):
                edge = (element_class.node1,element_class.node2)
                temp_graph = nx.contracted_edge( self.graph, edge, False )
        '''
    
    def find_fundamental_cycles(self):
        tmp_graph = nx.Graph(self.min_span_tree.copy()) # not a multigraph
        
        # iterating over each edge that does not belong to the spanning tree
        for element_class in self.flux_elements.values():
            # adding this edge to the spanning tree
            tmp_graph.add_edge(element_class.node1,element_class.node2,element=element_class)
            
            # finding fundamental cycle corresponding to the current edge
            paths = list(nx.all_simple_paths(tmp_graph,source=element_class.node1,target=element_class.node2))
            simple_cycle = paths[1] + list(reversed(paths[0][:-1]))
            # fundamental cycle orientation corresponds to the 
            # element_class.node1 -> element_class.node2 orientation of the cycle
            
            # generating list of elements in this             
            cycle_elements = []
            for i in range(len(simple_cycle)-1):
                element = tmp_graph[simple_cycle[i]][simple_cycle[i+1]]["element"]
                cycle_elements.append(element)
                
            self.fundamental_cycles_elements[element_class.refDes] = cycle_elements
            tmp_graph.remove_edge(element_class.node1,element_class.node2)

    def fundamental_cycles_info_str(self):        
        # iterating over fundamental cycles
        info_string = ""
        for flux_element_key,cycle_elements in self.fundamental_cycles_elements.items():
            element = self.elements[flux_element_key]
            info_string += "cycle for " + str(flux_element_key) + ", oriented as " + \
            str(element.node1) + "->" + str(element.node2) + ":\n"
            
            for element in cycle_elements:
                info_string += element.refDes + " "
                
            info_string.strip()
            info_string += '\n'
            
        return info_string
    
    def variables_info_str(self):
        info_str = ""
        for var in self.params.values():
            info_str += str(var.sym) + " = " + str(var.val) + "\n"
        return info_str
    
    def refresh_symbols(self):
        self.Hc_sym_phase = None
        self.Hj_sym_phase = None
        self.H_coupling_sym_phase = None
        self.H_external_sym_phase = None
        self._construct_Hc_sym_phase()
        self._construct_Hj_sym_phase()
        self._construct_H_coupling_sym_phase()
        self._construct_H_external_sym_phase()
        
        self.Hc_sym_cooperN = None
        self.Hj_sym_cooperN = None
        self.H_coupling_sym_cooperN = None
        self.H_external_sym_cooperN = None
        self._construct_Hc_sym_cooperN()
        self._construct_Hj_sym_cooperN()
        self._construct_H_coupling_sym_cooperN()
        self._construct_H_external_sym_cooperN()
        
        self.H_cooperN = self.Hc_sym_cooperN + self.Hj_sym_cooperN
                
        self.C_matrix_sym = None
        self._construct_Cap_matrix_sym_from_graph()
        
        self.C_matrix_intermidiate_sym = None
        self._construct_Cap_matrix_intermidiate_sym_from_Cap_matrix_sym()
        
        self.C_matrix_internal_sym = None
        self._construct_Cap_matrix_internal_sym_from_Cap_matrix_intermidiate_sym()
    
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
            
        # supply Elements with reference to this Scheme class
        self._provide_elements_with_scheme_instance()
        
        for element in self.elements.values():
            if( element.element_type == "EXT_V_SOURCE" ):
                self.ext_V_elements[element.refDes] = element
    
    
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
    
    
    def set_cooperN(self,cooper_N):
        self._construct_ops(cooper_N)
        self.cooper_N = cooper_N
    
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

    def _construct_H_coupling_sym_phase(self):
        pass

    def _construct_H_coupling_num_phase(self)             :
        raise NotImplementedError
        
    def _construct_H_coupling_sym_cooperN(self):
        pass 
    
    def _construct_H_coupling_num_cooperN(self, cooperN):
        raise NotImplementedError
        
    def _construct_H_external_sym_phase(self):
        pass

    def _construct_H_external_num_phase(self)             :
        raise NotImplementedError
        
    def _construct_H_external_sym_cooperN(self):
        pass
    
    def _construct_H_external_num_cooperN(self, cooperN):
        raise NotImplementedError
        
    
    
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
                
                # Zero row and column corresponds to the 
                # ground. This row and column is excluded
                if( node1 != 0 and node2 != 0 ):
                    C_matrix[node1-1,node2-1] -= C
                    C_matrix[node2-1,node1-1] -= C
                if( node1 != 0 ):
                    C_matrix[node1-1,node1-1] += C
                if( node2 != 0 ):
                    C_matrix[node2-1,node2-1] += C

        self.C_matrix_sym = C_matrix
    
    def _construct_Cap_matrix_intermidiate_sym_from_Cap_matrix_sym(self):
        '''
        @description:
            Constructs intermidiate matrix from full capcity matrix
        '''
                
        N = self.nodes_N - 1
        S_matrix = sympy.eye(N)        
        for ext_V_element in self.ext_V_elements.values():
            k = ext_V_element.node2 - 1
            l = ext_V_element.node1 - 1
            S_matrix += sympy.Matrix(N,N,lambda i,j: 1 if( i == k and j == l ) else 0)
            
        self.C_intermediate_sym = S_matrix.T * self.C_matrix_sym * S_matrix
        
    
    def _construct_Cap_matrix_internal_sym_from_Cap_matrix_intermidiate_sym(self):
        '''
        @description:
            Constructs internal 
        '''
            
        
        
                
    
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