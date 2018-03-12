# -*- coding: utf-8 -*-
from collections import OrderedDict
import itertools

import numpy as np

if( __name__ == "__main__" ):
    dict1 = OrderedDict( {"a":1.0,"b":2.0} )
    dict2 = OrderedDict( {"t":np.linspace(0.4,0.6,11)} )
    
    print( "dict1 items: ",list(dict1.values()) )
    print( "dict2 items: ",list(dict2.values()) )
    
    print( "product\n" )
    
    mesh = itertools.product( *list(dict2.keys()) ) 
    for point in mesh:
        print(point)
