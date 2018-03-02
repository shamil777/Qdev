from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from collections import OrderedDict

from .SLBase import SLBase

class ParametersSetupWidget(QtWidgets.QWidget,SLBase):
    def __init__(self, parent=None, flags=Qt.WindowFlags()):
        super(ParametersSetupWidget,self).__init__(parent,flags)
        SLBase.__init__(self)
        print("printing smth from ParametersSetupWidge")
        self.param_intervals_lineEdits = OrderedDict()
        self.param_intervals = OrderedDict()
        self.param_sym_labels = OrderedDict()
        self.colum_header_labels = []
        
        self.ref_to_parent = self.parent()
        
        self.init_GUI()
        
    def init_GUI(self):
        grid_layout = self.layout()

        if( grid_layout is None ):
            grid_layout = QtWidgets.QGridLayout()
            grid_layout.setSpacing(10)
            self.setLayout(grid_layout)
        
        params = self.ref_to_parent.simulator.scheme.params
        
        self.start_label = QtWidgets.QLabel(self)
        self.start_label.setText("start")
        self.stop_label = QtWidgets.QLabel(self)
        self.stop_label.setText("end")
        self.nPts_label = QtWidgets.QLabel(self)
        self.nPts_label.setText("N points")
        
        self.colum_header_labels = [self.start_label,self.stop_label,self.nPts_label]
        
        grid_layout.addWidget(self.start_label,0,1)
        grid_layout.addWidget(self.stop_label,0,2)
        grid_layout.addWidget(self.nPts_label,0,3)
        
        self.param_intervals = OrderedDict()
        for i,param in enumerate(params.values()):
            param_sym_label = QtWidgets.QLabel(self)
            param_sym_label.setText(str(param.sym))
            self.param_sym_labels[param.sym] = param_sym_label
            
            param_start_lineEdit = QtWidgets.QLineEdit(self)
            param_stop_lineEdit = QtWidgets.QLineEdit(self)
            param_nPts_lineEdit = QtWidgets.QLineEdit(self)
            
            param_start_lineEdit.editingFinished.connect(self.editing_finished_handler)
            param_stop_lineEdit.editingFinished.connect(self.editing_finished_handler)
            param_nPts_lineEdit.editingFinished.connect(self.editing_finished_handler)
            
            self.param_intervals_lineEdits[param.sym] = [param_start_lineEdit,param_stop_lineEdit,param_nPts_lineEdit]
            
            format_str = "{:.2g}"
            if( param.val is not None ):
                param_start_lineEdit.setText(format_str.format(param.val))
                param_stop_lineEdit.setText(format_str.format(param.val))
                param_nPts_lineEdit.setText(str(1))
                self.param_intervals[param.sym] = [param.val,param.val,1]
            else:
                self.param_intervals[param.sym] = [None]*3
                
            grid_layout.addWidget(param_sym_label,i+1,0)
            grid_layout.addWidget(param_start_lineEdit,i+1,1)
            grid_layout.addWidget(param_stop_lineEdit,i+1,2)
            grid_layout.addWidget(param_nPts_lineEdit,i+1,3)
        self.show()
    
    def refresh_data(self):
        print("refreshing data")
        
        layout = self.layout()

        for i in reversed(range(layout.count())): 
            layout.itemAt(i).widget().setParent(None)
            
        self.init_GUI()
            
    
    def editing_finished_handler(self):
        signal_source = self.sender()
        for i,sym in enumerate(self.param_intervals.keys()):
            for j in range(0,3):
                print("j = ",j)
                if( self.param_intervals_lineEdits[sym][j] == signal_source ):
                    if( j == 0 or j == 1 ):
                        self.param_intervals[sym][j] = float(signal_source.text())
                    else:
                        self.param_intervals[sym][j] = int(signal_source.text())
                        
                    return
                
    def transfer_internal_to_widget(self):
        print("transfer invoked in ParametersSetupWidget")
        