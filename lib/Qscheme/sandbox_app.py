import os
import sys

from Qscheme.simulator.SchemeSimulator import SchemeSimulator

if(__name__ == "__main__"):
    ss = SchemeSimulator()
    sys.exit(ss.GUI_loop())