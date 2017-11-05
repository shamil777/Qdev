# -*- coding: utf-8 -*-

import sympy

class Var():
    def __init__( self, sym=None,val=None ):
        if( not isinstance(sym,sympy.Symbol) ):
            sym = sympy.Symbol(sym)
        self.sym = sym
        self.val = val
        
def vars_from_symbols( symbols ):
    return [Var(symbol,None) for symbol in symbols]