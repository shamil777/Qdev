import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

from qscheme.variables import Var
import qscheme.SI as SI

import sympy
import numpy as np


class Element():
    def __init__(self, element_type, refDes, group_name, edge ):
        if( group_name == None ):
            group_name = refDes
            
        self.element_type = element_type # string that represents the element type
        self.group_name = group_name # user defined group name of the element
        self.refDes = refDes # unique designator of the element in scheme
        
        # whether or not this edge is not in the spanning tree
        self.flux_enabled = False
        
        self.scheme = None # reference to the scheme class
        
        # nodes that are the ends of the elements edge
        self.node1 = edge[0]
        self.node2 = edge[1]
        
        # common variables for every element
        self.f_ext = Var("\Phi_" + str(self.refDes) + "")        

        # list of Vars from variables.py of this element parameters
        self.params = [] 
        
        # subscript/designator of this element
        # for symbolic representation in formulas
        self.subscript = self.refDes
        self.f_ext_subscript = self.refDes
       
    # enables symbolic treatment of the external flux
    def flux_enable(self):
        self.flux_enabled = True
        self.params.append(self.f_ext)
        self.scheme.params[self.f_ext.sym] = self.f_ext
    
    def flux_disable(self):
        self.flux_enabled = False
        flux_idx = self.params.index(self.f_ext)
        self.params.pop(flux_idx)
        del self.scheme.params[self.f_ext]
        
        
    # changes subscript of the element and updates
    # all variables which symbolic representation is
    # dependent on this subscript
    def change_subscript(self,new_subscript, new_flux_subscript=None ):
        if( new_flux_subscript != None ):
            self.f_ext_subscript = new_flux_subscript
        self.subscript = new_subscript
        self._update_symbols()
    
    # updates symbols depending on the changes in
    # child-class data that affects symbolic 
    # representation
    def _update_symbols(self):
        raise NotImplementedError
    
    
    # This function must set class parameters
    # to refer to the Element.params list elements.
    # This function is used by Scheme class
    # right after collecting all variables of the scheme 
    # with unique names into Scheme.params list.
    # This function is used to force different Element class variables
    # with the same name
    # refer to the same variable structure, that can be accessed
    # either through Scheme.params, Element.params or Element."variable_name"
    # all of 3 ways to get access to the variables
    # become equal after this operation end in Scheme class
    def _connect_to_params_list(self):
        raise NotImplementedError
    
    # returns symbolic representations of kinetic and potential
    # hamiltonians in cooper number and phase basises
    def _Hc_symbol_cooperN(self):
        raise NotImplementedError
        
    def _Hc_symbol_phase(self):
        raise NotImplementedError
    
    def _Hj_symbol_cooperN(self):
        raise NotImplementedError
        
    def _Hj_symbol_phase(self):
        raise NotImplementedError
        
        
    # returns numerical representations of kinetic and potential
    # hamiltonians in cooper number and phase basises
    # objects returned as a QuTip operator classes
    def _Hc_num_cooperN(self, *ops):
        raise NotImplementedError
        
    def _Hc_num_phase(self, *ops):
        raise NotImplementedError
    
    def _Hj_num_cooperN(self, *ops):
        raise NotImplementedError
        
    def _Hj_num_phase(self, *ops):
        raise NotImplementedError
        
    
    # checking, whether the two Elements are equal
    # comparing they refDes members
    def __eq__( self, Element2 ):
        if( self.RefDes == Element2.RefDes ):
            return True
        else:
            return False
    
    # checking whether two Elements are NOT equal
    def __nq__( self, Element2 ):
        return (not self.__eq__( Element2 ))

# superconductor-insulator-superconductor josephson junction class
class SIS_JJ( Element ):
    def __init__(self, refDes, group_name,  edge ):     
        super( SIS_JJ,self ).__init__( "SIS_JJ", refDes, group_name, edge )
        # params with changable symbols
        self.C = Var("C_" + str(self.subscript) + "")
        self.Ej = Var("E_" + str(self.subscript) + "")
        
        # grouping params into list
        if( self.flux_enabled ):
            self.params = [self.C,self.Ej,self.f_ext]            
        else:
            self.params = [self.C, self.Ej]
    
    
    def _update_symbols(self):
        self.C.sym = sympy.Symbol("C_" + str(self.subscript) + "")
        self.Ej.sym = sympy.Symbol("E_" + str(self.subscript) + "")        
        if( self.flux_enabled ):
            self.f_ext.sym = sympy.Symbol("\Phi_" + str(self.f_ext_subscript) + "")


    def _connect_to_params_list(self):
        self.C = self.params[0]
        self.Ej = self.params[1]
        
        if( self.flux_enabled ):
            self.f_ext = self.params[2]
        
    
    def _Hc_symbol_phase(self):
        phi1_d = self.scheme.phi_d[self.node1].sym
        phi2_d = self.scheme.phi_d[self.node2].sym
        C = self.C.sym
            
        return C*SI.hbar.sym**2/(4*SI.e.sym**2)*(phi1_d - phi2_d)**2/2
    
    def _Hc_symbol_cooperN(self):
        n1 = self.scheme.n[self.node1].sym
        n2 = self.scheme.n[self.node2].sym
        C = self.C.sym
            
        return (n1-n2)**2*(4*SI.e.sym**2)/(2*C)
    
    def _Hj_symbol_phase(self):
        phi1 = self.scheme.phi[self.node1].sym
        phi2 = self.scheme.phi[self.node2].sym
        E_J = self.Ej.sym
        phi_ext = self.f_ext.sym
            
        if( self.flux_enabled is True ):
            return -E_J*sympy.cos(phi1 - phi2 + 2*sympy.pi*phi_ext)
        else:
            return -E_J*sympy.cos(phi1 - phi2)
        
    def _Hj_symbol_cooperN(self):
        l1 = self.scheme.l_ops[self.node1].sym
        l2 = self.scheme.l_ops[self.node2].sym
        r1 = self.scheme.r_ops[self.node1].sym
        r2 = self.scheme.r_ops[self.node2].sym
        
        E_j = self.Ej.sym
        phi_ext = self.f_ext.sym
        exp = sympy.exp(2*sympy.pi*sympy.I*phi_ext)
        if( self.flux_enabled is True ):
            return -E_j*(l1*r2*exp + r1*l2/exp)/2    
        else:
            return -E_j*(l1*r2 + r1*l2)/2
        
    def _Hj_num_cooperN(self):      
        Ej = self.Ej.val
        Hj_num_cooperN = 0
        
        exp = 1
        if( self.flux_enabled == True ):
            exp = np.exp(1j*2*np.pi*self.f_ext.val)
        
        node1 = self.node1
        node2 = self.node2
        l1 = self.scheme.l_ops[node1].val
        l2 = self.scheme.l_ops[node2].val
        r1 = self.scheme.r_ops[node1].val
        r2 = self.scheme.r_ops[node2].val
        
        Hj_num_cooperN += -Ej/2*(r1*l2*exp + l1*r2/exp)
        return Hj_num_cooperN

        
        
        
class Cap( Element ):
    def __init__(self, refDes, group_name, edge, capacitance=51e-15 ):                                        
        super( Cap,self ).__init__( "Cap", refDes, group_name, edge )
        
        self.subscript = str(self.node1) + str(self.node2)
        self.C = Var("C_" + str(self.subscript) + "")
        self.params = [self.C]

            
    
    def _update_symbols(self):
        self.C.sym = sympy.Symbol("C_" + str(self.subscript) + "")
        
    def _connect_to_params_list(self):
        self.C = self.params[0]
    
    def _Hc_symbol_phase(self):
        phi1_d = self.scheme.phi_d[self.node1].sym
        phi2_d = self.scheme.phi_d[self.node2].sym
        
        C = self.C.sym
            
        return C*SI.hbar.sym**2/(4*SI.e.sym**2)*(phi1_d - phi2_d)**2/2
    
    def _Hc_symbol_cooperN(self):
        n1 = self.scheme.n[self.node1].sym
        n2 = self.scheme.n[self.node2].sym

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

        
class ExternalForceElement( Element ):
    def __init__(self, element_type, refDes,group_name,edge ):
        super( ExternalForceElement,self ).__init__( element_type, refDes, group_name, edge )

class Battery( ExternalForceElement ):
    def __init__(self, refDes, group_name, edge ):
        super( Battery,self ).__init__( "EXT_V_SOURCE", refDes, group_name, edge )
    
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
    
class Resonator_V( ExternalForceElement ):
    def __init__(self, refDes, group_name, edge ):
        super( Resonator_V,self ).__init__( "EXT_V_SOURCE", refDes, group_name, edge )
    
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
    