# -*- coding: utf-8 -*-
from collections import OrderedDict
import itertools

import numpy as np

if( __name__ == "__main__" ):
    dict1 = OrderedDict( {"a":1.0,"b":2.0} )
    dict2 = OrderedDict( {"t":np.linspace(0.4,0.6,3),"t2":np.linspace(10,20,3)} )
    for key,val in dict1.items():
        dict1[key] = [val]
        
    premesh_sweep = itertools.product(*[itertools.product(key,val) for key,val in dict(**dict1,**dict2).items()])
    
    for point in premesh_sweep:
        print(point)
    