import csv
import re
import networkx as nx

import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

from Qscheme.variables import Var,vars_from_symbols
import Qscheme.SI as SI

import sympy
import numpy as np
import qutip as qp
from sympy import init_printing
init_printing()

from collections import OrderedDict


class Element():
    def __init__(self, element_type, refDes, group_name, edge ):
        if( group_name == None ):
            group_name = refDes
        self.element_type = element_type
        self.group_name = group_name
        self.refDes = refDes
        
        self.node1 = edge[0]
        self.node2 = edge[1]

        self.params = []      
        self.subscript = None
       
    def change_subscript(self,new_subscript):
        self.subscript = new_subscript
        self._update_symbols()
    
    def _init_symbols(self):
        raise NotImplementedError
    
    def _update_symbols(self):
        raise NotImplementedError
        
    def _connect_to_params_list(self):
        raise NotImplementedError
    
    def _Hc_symbol_cooperN(self):
        raise NotImplementedError
        
    def _Hc_symbol_phase(self):
        raise NotImplementedError
    
    def _Hj_symbol_cooperN(self):
        raise NotImplementedError
        
    def _Hj_symbol_phase(self):
        raise NotImplementedError
        
    def _Hc_num_cooperN(self, *ops):
        raise NotImplementedError
        
    def _Hc_num_phase(self, *ops):
        raise NotImplementedError
    
    def _Hj_num_cooperN(self, *ops):
        raise NotImplementedError
        
    def _Hj_num_phase(self, *ops):
        raise NotImplementedError
        
    
    
    def __eq__( self, Element2 ):
        if( self.RefDes == Element2.RefDes ):
            return True
        else:
            return False
        
    def __nq__( self, Element2 ):
        return (not self.__eq__( Element2 ))
        
class SIS_JJ( Element ):
    def __init__(self, refDes, group_name,  edge ):     
        super( SIS_JJ,self ).__init__( "SIS_JJ", refDes, group_name, edge )
        self.subscript = refDes
        
        # params with changable symbols
        self.C = Var("C_{" + str(self.subscript) + "}")
        self.Ej = Var("E_{" + str(self.subscript) + "}")
        self.f_ext = Var("\\varphi_{" + str(self.refDes) + "}")        
        self.params = [self.C,self.Ej,self.f_ext]
        
        # params with fixed symbols 
        self.n1,self.n2 = vars_from_symbols( 
                        sympy.symbols(      
                                        ["\hat{n}_{"+str(self.node1) + "}",
                                         "\hat{n}_{"+str(self.node2) + "}"],
                                        commutative=False
                                    ) 
                        )
        self.phi1_d,self.phi2_d = vars_from_symbols(
                            sympy.symbols(
                                        ["\dot{\hat{\\varphi}}_" + str(self.node1),
                                         "\dot{\hat{\\varphi}}_" + str(self.node2)],
                                         commutative=False
                                         )
                            )
        self.phi1, self.phi2 = vars_from_symbols(
                            sympy.symbols(
                                        ["\hat{\\varphi}_{" + str(self.node1) + "}",
                                         "\hat{\\varphi}_{" + str(self.node2) + "}"]
                                        )
                            )
        
        self.l1,self.l2,self.r1,self.r2 = vars_from_symbols(
                            sympy.symbols(
                                        ["\hat{b}_{"+str(self.node1) + "}",
                                         "\hat{b}_{"+str(self.node2)+"}",
                                         "\hat{b}^\dag_{"+str(self.node1) + "}",
                                         "\hat{b}^\dag_{"+str(self.node2) + "}"],
                                        commutative = False
                                        )
                            )        
        
        self.ext_flux_subscript = refDes
        
        self.flux_enabled = False
    
    def _update_symbols(self):
        self.C.sym = sympy.Symbol("C_{" + str(self.subscript) + "}")
        self.Ej.sym = sympy.Symbol("E_{" + str(self.subscript) + "}")        
        self.f_ext.sym = sympy.Symbol("\\varphi_{" + str(self.ext_flux_subscript) + "}")
        
    def change_subscript(self,new_subscript,new_flux_subscript=None):
        if( new_flux_subscript != None ):
            self.ext_flux_subscript = new_flux_subscript
        super(SIS_JJ,self).change_subscript(new_subscript)
        
    def _connect_to_params_list(self):
        self.C = self.params[0]
        self.Ej = self.params[1]
        self.f_ext = self.params[2]
        
    def enable_ext_flux(self):
        self.flux_enabled = True
    
    def _Hc_symbol_phase(self):
        phi1_d,phi2_d = self.phi1_d.sym,self.phi2_d.sym
        C = self.C.sym
            
        return C*SI.hbar.sym**2/(4*SI.e.sym**2)*(phi1_d - phi2_d)**2/2
    
    def _Hc_symbol_cooperN(self):
        n1 = self.n1.sym
        n2 = self.n2.sym
        C = self.C.sym
            
        return (n1-n2)**2*(4*SI.e.sym**2)/(2*C)
    
    def _Hj_symbol_phase(self):
        phi1,phi2 = self.phi1.sym, self.phi2.sym
        E_J = self.Ej.sym
        phi_ext = self.f_ext.sym
            
        if( self.flux_enabled is True ):
            return -E_J*sympy.cos(phi1 - phi2 + 2*sympy.pi*phi_ext)
        else:
            return -E_J*sympy.cos(phi1 - phi2)
        
    def _Hj_symbol_cooperN(self):
        l1,l2,r1,r2 = self.l1.sym,self.l2.sym,self.r1.sym,self.r2.sym
        E_j = self.Ej.sym
        phi_ext = self.f_ext.sym
        
        if( self.flux_enabled is True ):
            return -E_j*(l1*r2*sympy.exp(1j*2*sympy.pi*phi_ext) + r1*l2*sympy.exp(-1j*2*sympy.pi*phi_ext))/2    
        else:
            return -E_j*(l1*r2 + r1*l2)/2
        
    def _Hj_num_cooperN(self,raising_ops,lowering_ops):        
        node1 = self.node1
        node2 = self.node2
        Ej = self.Ej.val
        Hj_num_cooperN = 0*raising_ops[0]
        
        exp = 1
        if( self.flux_enabled == True ):
            exp = np.exp(1j*2*np.pi*self.f_ext.val)
        
        if( self.node1 != 0 and self.node2 != 0 ):
            Hj_num_cooperN += -Ej/2*(raising_ops[node1-1]*lowering_ops[node2-1]*exp \
                                          + lowering_ops[node1-1]*raising_ops[node2-1]/exp)
        else:
            node = node1 + node2 # equals to nonzero node number
            Hj_num_cooperN += -Ej/2*(lowering_ops[node-1]*exp \
                                          + raising_ops[node-1]/exp)
            
        return Hj_num_cooperN

        
        
        
class Cap( Element ):
    def __init__(self, refDes, group_name, edge, capacitance=51e-15 ):                                        
        super( Cap,self ).__init__( "Cap", refDes, group_name, edge )
        
        self.subscript = str(self.node1) + str(self.node2)
        self.C = Var("C_{" + str(self.subscript) + "}")
        self.params = [self.C]
        
        self.n1,self.n2 = vars_from_symbols( 
                        sympy.symbols(      
                                        ["\hat{n}_{"+str(self.node1) + "}",
                                         "\hat{n}_{"+str(self.node2) + "}"],
                                        commutative=False
                                    ) 
                        )
        self.phi1_d,self.phi2_d = vars_from_symbols(
                            sympy.symbols(
                                        ["\dot{\hat{\\varphi}}_" + str(self.node1),
                                         "\dot{\hat{\\varphi}}_" + str(self.node2)],
                                         commutative=False
                                         )
                            )

            
    
    def _update_symbols(self):
        self.C.sym = sympy.Symbol("C_{" + str(self.subscript) + "}")
        
    def _connect_to_params_list(self):
        self.C = self.params[0]
    
    def _Hc_symbol_phase(self):
        phi1_d,phi2_d = self.phi1_d.sym,self.phi2_d.sym
        
        C = self.C.sym
            
        return C*SI.hbar.sym**2/(4*SI.e.sym**2)*(phi1_d - phi2_d)**2/2
    
    def _Hc_symbol_cooperN(self):
        n1,n2 = self.n1.sym,self.n2.sym

        C = self.C.sym
            
        return (n1-n2)**2*(4*SI.e.sym**2)/(2*C)
    
    def _Hj_symbol_phase(self):
        return 0
    
    def _Hj_symbol_cooperN(self):
        return 0
    
    def _Hc_num_cooperN(self, *ops):
        return 0
        
    def _Hc_num_phase(self, *ops):
        return 0
    
    def _Hj_num_cooperN(self, *ops):
        return 0
        
    def _Hj_num_phase(self, *ops):
        return 0

        

class Battery( Element ):
    def __init__(self, refDes, group_name, edge ):
        super( Battery,self ).__init__( "Battery", refDes, group_name, edge )
    
    def _init_symbols(self):
        pass
    
    def _update_symbols(self):
        pass
    
    def _connect_to_params_list(self):
        pass
    
    def _Hc_symbol_phase(self):
        return 0
    
    def _Hc_symbol_cooperN(self):
        return 0
    
    def _Hj_symbol_phase(self):
        return 0
    
    def _Hj_symbol_cooperN(self):
        return 0

    def _Hc_num_cooperN(self, *ops):
        return 0
        
    def _Hc_num_phase(self, *ops):
        return 0
    
    def _Hj_num_cooperN(self, *ops):
        return 0
        
    def _Hj_num_phase(self, *ops):
        return 0

def _get_refDes_groupNames_PADS( file_rows ):
    ''' 
    @description:
        parses file and return edge types, represented by electrical elements
    @parameters:
        file_rows - list of lists containing csv file (e.g. list( csv.reader(file)) 
    @return: list of Element() objects used in this netlist file
    '''
    refDes_list = []
    group_name_list = []
    
    element_i = -1
    for i,row in enumerate( file_rows ):
        if( row[0] == "*PART*" ):
            element_i = i+1
            break
    
    row = file_rows[element_i]
    while( len(row) > 0 ):
        refDes = row[0]
        group_name = row[5].split(',')[1]
        
        refDes_list.append(refDes)
        group_name_list.append(group_name)
        
        
        element_i += 1
        row = file_rows[element_i]
        
    return refDes_list, group_name_list
	
def build_graph_from_netlist_PADS( file_rows ):
    ''' 
    @description:
        parses file and MultiGraph object that represents the schematic
        Nodes are nets in electrical schematic, enumeration is kept the same as in the file
        Edges are electrical elements (capacitors, SIS josephson junctions etc.)
    @parameters:
        file_rows - list of lists containing csv file (e.g. list( csv.reader(file)) 
    @return: MultiGraph object from networkx module, containing electrical schematic
    '''
    refDes_list, group_name_list = _get_refDes_groupNames_PADS( file_rows )        
    nets_list = [[] for x in refDes_list]
    graph = nx.MultiGraph()
    
    for i,row in enumerate( file_rows ):
        if( len(row) < 2 ):
            continue
            
        net_id = row[1].split('_') 
        
        # found row with net information (net is the node of the graph)
        if( net_id[0] == "Net" ):
            graph.add_node( int(net_id[1]) ) # adding net with corresponding number
            element_row = file_rows[i+1]
            for element in element_row[:-1]: # last element is empty due to DipTrace export...
                refDes,pinDesc = element.split('.')
                i = refDes_list.index( refDes )
                if( pinDesc == "NEG" ):
                    nets_list[i].insert( 0,int(net_id[1]) )
                else:
                    nets_list[i].append( int(net_id[1]) )
    
    for i,edge in enumerate(nets_list):
        refDes = refDes_list[i]
        element_arg = None
        refName = re.compile("[A-Za-z]*").match(refDes)[0]
        group_name = group_name_list[i]
        
        if( refName == "J" ): # Josephson junction
            element_arg = SIS_JJ(refDes,group_name,edge)
        elif( refName == "C" ): #capacitor
            element_arg = Cap(refDes,group_name, edge)
        elif( refName == "B" ): # battery
            element_arg = Battery(refDes,group_name,edge)
          
        graph.add_edge( *edge, element=element_arg )
        
    return graph

def build_graph_from_file_PADS( file_path ):
    with open( file_path, "r" ) as file:
        rows = list(csv.reader(file, delimiter=' '))
    return build_graph_from_netlist_PADS( rows )