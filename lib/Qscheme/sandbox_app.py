import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

from Qscheme.simDialog import SimWindow

from PyQt5 import QtWidgets

from sympy import init_printing
init_printing()

if( __name__ == "__main__" ):
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()     

    window = SimWindow()
    sys.exit(app.exec_())