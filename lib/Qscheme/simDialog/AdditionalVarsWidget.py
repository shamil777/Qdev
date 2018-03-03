from PyQt5 import QtCore,QtWidgets
from PyQt5.QtCore import Qt

from .SLBase import SLBase

class AdditionalVarsWidget(QtWidgets.QWidget,SLBase):
    def __init__(self, parent=None, flags=Qt.WindowFlags()):
        super(AdditionalVarsWidget,self).__init__(parent,flags)      
        
        self._fill_SL_dicts()
        self.init_GUI()
        
    def init_GUI(self):
        grid = QtWidgets.QGridLayout()
        grid.setSpacing(10)
        self.setLayout(grid)
        
        for i in range(0,5):
            symbol_str = QtWidgets.QLineEdit(self)
            grid.addWidget(symbol_str,i,0)
            val_str = QtWidgets.QLineEdit(self)
            grid.addWidget(val_str,i,1)
            
        self.show()
    
    def _fill_SL_dicts(self):
        pass
    
    def transfer_internal_to_widget(self):
        print("transfer invoked in AdditionalVarsWidget")
        
        