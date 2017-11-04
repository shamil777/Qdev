import csv
import re
import networkx as nx

import sympy
import numpy as np
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

        self.params = OrderedDict()        
        self.subscript = None
        self._init_symbols()
       
    def change_subscript(self,new_subscript):
        self.subscript = new_subscript
        self._update_symbols()
    
    def _init_symbols(self):
        raise NotImplementedError
    
    def _update_symbols(self):
        raise NotImplementedError
    
    def _Hc_symbol_cooperN(self):
        raise NotImplementedError
        
    def _Hc_symbol_phase(self):
        raise NotImplementedError
    
    def _Hj_symbol_cooperN(self):
        raise NotImplementedError
        
    def _Hj_symbol_phase(self):
        raise NotImplementedError
    
    def __eq__( self, Element2 ):
        if( self.RefDes == Element2.RefDes ):
            return True
        else:
            return False
        
    def __nq__( self, Element2 ):
        return (not self.__eq__( Element2 ))
        
class SIS_JJ( Element ):
    def __init__(self, refDes, group_name,  edge, capacitance=4.5e-15, E_J = 38e9 ):     
        self.C = None
        self.Ej = None
        self.f_ext = None
        self.ext_flux_subscript = None
        
        self.flux_enabled = False
        
        super( SIS_JJ,self ).__init__( "SIS_JJ", refDes, group_name, edge )
    
    def _init_symbols(self):
        if(self.subscript is None):
            self.subscript = self.refDes
        self.C = sympy.Symbol("C_{" + str(self.subscript) + "}")
        self.Ej = sympy.Symbol("E_{" + str(self.subscript) + "}")        
        self.f_ext = sympy.Symbol("\\varphi_{ext}")
        self.params[self.C] = None
        self.params[self.Ej] = None
        self.params[self.f_ext] = None
    
    def _update_symbols(self):
        Cj_val = self.params[self.C]
        Ej_val = self.params[self.Ej]
        f_ext_val = self.params[self.f_ext]
        self.C = sympy.Symbol("C_{" + str(self.subscript) + "}")
        self.Ej = sympy.Symbol("E_{" + str(self.subscript) + "}")        
        self.f_ext = sympy.Symbol("\\varphi_{" + str(self.ext_flux_subscript) + "}")
        self.params = OrderedDict()
        self.params[self.C] = Cj_val
        self.params[self.Ej] = Ej_val
        self.params[self.f_ext] = f_ext_val
        
    def change_subscript(self,new_subscript,new_flux_subscript="ext"):
        self.ext_flux_subscript = new_flux_subscript
        super(SIS_JJ,self).change_subscript(new_subscript)
        
    def enable_ext_flux(self):
        self.flux_enabled = True
    
    def _Hc_symbol_phase(self):
        phi1_d,phi2_d,hbar,e = sympy.symbols(("\dot{\hat{\\varphi}}_" + str(self.node1),
                                             "\dot{\hat{\\varphi}}_" + str(self.node2),
                                             "\hbar","e") )
        C = self.C
            
        return C*hbar**2/(4*e**2)*(phi1_d - phi2_d)**2/2
    
    def _Hc_symbol_cooperN(self):
        n1,n2,e = sympy.symbols(["\hat{n}_{"+str(self.node1) + "}",
                               "\hat{n}_{"+str(self.node2) + "}",
                               "e"])
        C = self.C
            
        return (n1-n2)**2*(4*e**2)/(2*C)
    
    def _Hj_symbol_phase(self):
        phi1, phi2 = sympy.symbols(["\hat{\\varphi}_{" + str(self.node1) + "}",
                                    "\hat{\\varphi}_{" + str(self.node2) + "}"])
        E_J = self.Ej
        phi_ext = self.f_ext
            
        if( self.flux_enabled is True ):
            return -E_J*sympy.cos(phi1 - phi2 + 2*sympy.pi*phi_ext)
        else:
            return -E_J*sympy.cos(phi1 - phi2)
        
    def _Hj_symbol_cooperN(self):
        l1,l2,r1,r2 = sympy.symbols( ["\hat{b}_{"+str(self.node1) + "}",
                                 "\hat{b}_{"+str(self.node2)+"}",
                                 "\hat{b}^\dag_{"+str(self.node1) + "}",
                                 "\hat{b}^\dag_{"+str(self.node2) + "}"])
        E_j = self.Ej
        phi_ext = self.f_ext
        
        if( self.flux_enabled is True ):
            return -E_j*(l1*r2*sympy.exp(1j*2*sympy.pi*phi_ext) + r1*l2*sympy.exp(-1j*2*sympy.pi*phi_ext))/2    
        else:
            return -E_j*(l1*r2 + r1*l2)/2
        
        
        
class Cap( Element ):
    def __init__(self, refDes, group_name, edge, capacitance=51e-15 ):                                        
        self.C = None
        
        super( Cap,self ).__init__( "Cap", refDes, group_name, edge )
    
    def _init_symbols(self):
        if( self.subscript is None ):
            self.subscript = str(self.node1) + str(self.node2)
        self.C = sympy.Symbol("C_{" + str(self.subscript) + "}")
        self.params[self.C] = None
    
    def _update_symbols(self):
        C_val = self.params[self.C]
        self.C = sympy.Symbol("C_{" + str(self.subscript) + "}")
        self.params[self.C] = C_val
    
    def _Hc_symbol_phase(self):
        phi1_d,phi2_d,hbar,e = sympy.symbols(("\dot{\hat{\\varphi}}_" + str(self.node1),
                                             "\dot{\hat{\\varphi}}_" + str(self.node2),
                                             "\hbar","e") )
        
        C = self.C
            
        return C*hbar**2/(4*e**2)*(phi1_d - phi2_d)**2/2
    
    def _Hc_symbol_cooperN(self):
        n1,n2,e = sympy.symbols(["\hat{n}_{"+str(self.node1) + "}",
                               "\hat{n}_{"+str(self.node2) + "}",
                               "e"])

        C = self.C
            
        return (n1-n2)**2*(4*e**2)/(2*C)
    
    def _Hj_symbol_phase(self):
        return 0
    
    def _Hj_symbol_cooperN(self):
        return 0
        

class Battery( Element ):
    def __init__(self, refDes, group_name, edge ):
        super( Battery,self ).__init__( "Battery", refDes, group_name, edge )
    
    def _init_symbols(self):
        pass
    
    def _update_symbols(self):
        pass
    
    def _Hc_symbol_phase(self):
        return 0
    
    def _Hc_symbol_cooperN(self):
        return 0
    
    def _Hj_symbol_phase(self):
        return 0
    
    def _Hj_symbol_cooperN(self):
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
                refDes = element.split('.')[0]
                i = refDes_list.index( refDes )
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