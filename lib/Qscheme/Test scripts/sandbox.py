from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

import sys

class myMain(QtWidgets.QMainWindow):
    def __init__(self):
        super(myMain,self).__init__()
        self.child_window = myChild(self)
        self.show()
        self.child_window.show()
        
class myChild(QtWidgets.QWidget):
    def __init__(self,parent=None):
        super(myChild,self).__init__(parent,Qt.Dialog)

if( __name__ == "__main__" ):
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()
    
    window = myMain()
    app.exec_()