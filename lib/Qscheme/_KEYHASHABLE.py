# -*- coding: utf-8 -*-
"""
Created on Tue May  8 14:28:02 2018

@author: SUNSHINE_MACHINE
"""

class SETTINGS_FIELD:
    LAST_QDEV_FILE_OPENED = "LAST_QDEV_FILE_OPENED"

class VAR_KIND:
    FIXED = "FIXED"
    SWEEP = "SWEEP"
    EQUATION = "EQUATION"
    
class SIM_SUBSYS:
    INTERNAL = "Internal"
    COUPLING = "Coupling"
    WHOLE = "Whole"
    
class SIM_BASIS:
    COOPER = "Cooper"
    PHASE = "Phase"
    
class KW:
    EVALS = "EVALS"
    EVECTS = "EVECTS"
    POINT = "POINT"
    
class IDX:
    POINT = 0
    EVALS = 1
    EVECTS = 2
    
class DRIVE_KWDS:
    pass    