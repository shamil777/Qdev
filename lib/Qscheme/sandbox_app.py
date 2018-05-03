import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\" )

from Qscheme.simDialog import SimWindow

from PyQt5.QtWidgets import QApplication
5
if( __name__ == "__main__" ):
    if not QApplication.instance():
        app = QApplication(sys.argv)
    else:
        app = QApplication.instance()     

    window = SimWindow()
    sys.exit(app.exec_())