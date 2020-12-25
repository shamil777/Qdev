# -*- coding: utf-8 -*-
import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )
from qscheme.variables import Var
import numpy as np

e = Var("e",1.6e-19) # electron charge SI
hbar = Var("\hbar",1.055e-34) # plank constant SI (this is \hbar)
h = Var("h",hbar.val*(2*np.pi))

alpha = 0.43 # ratio of the smaller and larger josephson junctions areas
E_C = 17.6 # GHz, roughly corresponding to C = 4.5 fF
C = 4*e.val**2/2/(h.val*E_C*10**9)
E_J = 36.2/alpha # GHz
Csh = 51e-15/C # in units of C (C equals to larger josephson junctions capacitance)

Z_r = 50 # Ohm - resonator line impendance
freq_r = 6.00 # GHz resonator frequency
omerga_r = freq_r*(2*np.pi)
beta = 0.043 # coupling parameter
#rms single-photon voltage in resonator
V_r0 = freq_r*np.sqrt(2*h.val*Z_r)