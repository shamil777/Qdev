from PyQt5 import QtCore,QtWidgets
from PyQt5.QtCore import Qt

from .SLBase import SLBase

class SimulationRegimeWidget(QtWidgets.QWidget,SLBase):
    def __init__(self, parent=None, flags=Qt.WindowFlags()):
        super(SimulationRegimeWidget,self).__init__(parent,flags)
        SLBase.__init__(self)
        
        self.cooperN = 5
        
        self.init_GUI()
        
    def init_GUI(self):
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.setSpacing(10)
        self.setLayout(v_layout)
        
        # 1st from top
        regime_label = QtWidgets.QLabel(self)
        regime_label.setText("Simulation context")
        
        # 2nd from top
        subsys_comboBox_layout = QtWidgets.QHBoxLayout()
        
        subsys_label = QtWidgets.QLabel("subsystem: ",self)
        subsys_comboBox = QtWidgets.QComboBox()
        subsys_comboBox.addItems(["Internal","External","Internal coupling","External coupling","All"])
        
        subsys_comboBox_layout.addWidget(subsys_label)
        subsys_comboBox_layout.addWidget(subsys_comboBox)
        
        # 3rd from top
        simType_hLayout = QtWidgets.QHBoxLayout()
        
        simType_label = QtWidgets.QLabel("type: ",self)
        simType_comboBox = QtWidgets.QComboBox(self)
        simType_comboBox.addItem("Single point")
        simType_comboBox.addItem("Parameters product sweep")
        
        simType_hLayout.addWidget(simType_label)
        simType_hLayout.addWidget(simType_comboBox)
        
        # 4th from top
        cooper_N_hLayout = QtWidgets.QHBoxLayout()
        
        cooper_N_label = QtWidgets.QLabel(self)
        cooper_N_label.setText("cooper pairs max")
                
        self.cooper_N_LineEdit = QtWidgets.QLineEdit(self)
        self.cooper_N_LineEdit.setText(str(self.cooperN))
        self.cooper_N_LineEdit.editingFinished.connect(self.cooper_N_LineEdit_editing_finished)
        
        cooper_N_hLayout.addWidget(cooper_N_label)
        cooper_N_hLayout.addWidget(self.cooper_N_LineEdit)
        
        # grouping all in v_laout
        v_layout.addWidget(regime_label)
        v_layout.setAlignment(regime_label,Qt.AlignTop)
        
        v_layout.addLayout(subsys_comboBox_layout)
        v_layout.setAlignment(subsys_comboBox_layout,Qt.AlignTop)
        
        v_layout.addLayout(simType_hLayout)
        v_layout.setAlignment(simType_hLayout,Qt.AlignTop)
        
        v_layout.addLayout(cooper_N_hLayout)
        v_layout.setAlignment(cooper_N_hLayout,Qt.AlignTop)
        self.show()
            
        
    def cooper_N_LineEdit_editing_finished(self):
        self.cooperN = int( self.cooper_N_LineEdit.text() )
        
    def transfer_internal_to_widget(self):
        print("transfer invoked in SimulationRegimeWidget")
            
        
        