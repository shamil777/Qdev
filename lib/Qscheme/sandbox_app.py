import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )


from Qscheme.simulator.SchemeSimulator import SchemeSimulator

if( __name__ == "__main__" ):
    ss = SchemeSimulator()
    sys.exit(ss.GUI_loop())
    