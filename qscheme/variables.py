# -*- coding: utf-8 -*-

import sympy

from .simulator import operator_constructors as oc

class Var():
    def __init__( self, sym=None,val=None ):
        if( not isinstance(sym,sympy.Symbol) ):
            try:
                sym = sympy.Symbol(sym)
            except TypeError as e:
                sym  = sympy.sympify(sym)
                
        self.sym = sym
        self.val = val
        
def vars_from_symbols( symbols ):
    return [Var(symbol,None) for symbol in symbols]

def _vars_to_dict( vars_list ):
    return {var.sym:var.val for var in vars_list}

def _get_cooperN_ops_as_vars(nodes_N,cooper_N=None):
    if( cooper_N is None ):
        charge_ops = [None for i in range(nodes_N)]
    else:
        charge_ops = oc.charge_ops(nodes_N - 1,cooper_N)
        
    charge_vars = [Var(0,0)]
    for i in range(nodes_N - 1):
        charge_vars.append( Var(sympy.Symbol("\hat{n}_{" + str(i+1) + "}",commutative=False),charge_ops[i]) )
        
    return charge_vars

def _get_raising_lowering_ops_as_vars(nodes_N, cooper_N=None):
    if( cooper_N is None ):
        raising = [None for i in range(nodes_N)]
        lowering = raising
    else:
        raising,lowering = oc.raising_lowering_ops(nodes_N-1, cooper_N)
        
    raising_vars = [Var(1,1)]
    lowering_vars = [Var(1,1)]
    for i in range(nodes_N - 1):
        raising_vars.append( Var(sympy.Symbol("\hat{b}^\dag_{" + str(i+1) + "}",commutative=False),raising[i]) )
        lowering_vars.append( Var(sympy.Symbol("\hat{b}_{" + str(i+1) + "}",commutative=False), lowering[i]) )
    return raising_vars,lowering_vars
    

def _get_phase_ops_as_vars(nodes_N):
    return [Var(sympy.Symbol("\hat{\varphi}_{" + str(i) + "}", commutative=False),None) for i in range(nodes_N)]

def _phase_pd_ops_as_vars(nodes_N):
    return [Var(sympy.Symbol("\hat{\dot{\varphi}}_{" + str(i) + "}", commutative=False),None) for i in range(nodes_N)]
