from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from collections import OrderedDict

from ..SLBase import SLBaseWidget

class SubscriptsWidget(QtWidgets.QWidget,SLBaseWidget):
    def __init__(self, parent=None, flags=Qt.WindowFlags()):
        super(SubscriptsWidget,self).__init__(parent,flags)
        
        self.ref_to_parent = self.parent()
        
        self.group_name_subscripts = OrderedDict()        
        self.subscript_lineEdits = OrderedDict()
        
        self.fill_SL_names()
        
    def init_GUI(self):
        v_layout = QtWidgets.QVBoxLayout()
        self.setLayout(v_layout)
        
        grid = QtWidgets.QGridLayout()
        grid.setSpacing(10)
        
        elements = self.ref_to_parent.simulator.scheme.elements
        
        for element in elements.values():
            self.group_name_subscripts[element.group_name] = ""
        
        for i,group_name in enumerate(self.group_name_subscripts.keys()):
            group_label = QtWidgets.QLabel(group_name,self)
            grid.addWidget(group_label,i,0)
            
            lineEdit = QtWidgets.QLineEdit(self)
            lineEdit.editingFinished.connect(self.editing_finished_handler)
            self.subscript_lineEdits[group_name] = lineEdit
            grid.addWidget(lineEdit,i,1)
        
        v_layout.addLayout(grid)
        
        assign_button = QtWidgets.QPushButton("Assign",self)
        assign_button.clicked.connect(self.assign_button_clicked_handler)
        
        v_layout.addWidget(assign_button)
        
        self.show()
        
    
    def fill_SL_names(self):
        self.SL_attributes_names = ["group_name_subscripts"]
    
    def assign_button_clicked_handler(self):
        scheme = self.ref_to_parent.simulator.scheme
        scheme.assign_subscripts_to_nameGroups(list(self.group_name_subscripts.keys()),list(self.group_name_subscripts.values()))
        self.ref_to_parent.parameter_setup_widget.scheme_subscripts_changed_handler()
        self.ref_to_parent.fit_window_to_content()
        
    def editing_finished_handler(self):        
        signal_source = self.sender()       
        for group_name in self.group_name_subscripts.keys():
            if( signal_source == self.subscript_lineEdits[group_name] ):
                self.group_name_subscripts[group_name] = signal_source.text()
                return
            
    def transfer_internal_to_widget(self):
        for key,value in self.subscript_lineEdits.items():
            value.setText(self.group_name_subscripts[key])
        
        