from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from ..SLBase import SLBaseWidget
from .._KEYHASHABLE import SIM_BASIS,SIM_SUBSYS

class SimulationRegimeWidget(QtWidgets.QWidget,SLBaseWidget):
    
    def __init__(self, parent=None, flags=Qt.WindowFlags()):
        super(SimulationRegimeWidget,self).__init__(parent,flags)
        
        # Save/Load variables section START #    
        self.cooperN = 5
        self.simulation_basis = None
        self.subsystem_choice = None
        # Save/Load variables section END #
        
        self.fill_SL_names()
        
        self.subsys_comboBox = QtWidgets.QComboBox()
        self.subsys_comboBox.addItems([SIM_SUBSYS.INTERNAL,
                                       SIM_SUBSYS.COUPLING,
                                       SIM_SUBSYS.WHOLE])
        self.subsys_comboBox.currentIndexChanged.connect(self.subsys_comboBox_changed_handler)
        self.subsys_comboBox_changed_handler()
        
        self.basis_comboBox = QtWidgets.QComboBox(self)
        self.basis_comboBox.addItems([SIM_BASIS.COOPER,
                                      SIM_BASIS.PHASE])
        self.basis_comboBox.currentIndexChanged.connect(self.basis_comboBox_changed_handler)
        self.basis_comboBox_changed_handler()
        
        self.cooper_N_LineEdit = QtWidgets.QLineEdit(self)
        self.cooper_N_LineEdit.setText(str(self.cooperN))
        self.cooper_N_LineEdit.editingFinished.connect(self.cooper_N_LineEdit_editing_finished)
        
    def init_GUI(self):
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.setAlignment(Qt.AlignTop)
        v_layout.setSpacing(10)
        self.setLayout(v_layout)
        
        # 1st from top
        regime_label = QtWidgets.QLabel(self)
        regime_label.setText("Simulation context")
        
        # 2nd from top
        subsys_comboBox_layout = QtWidgets.QHBoxLayout()
        
        subsys_label = QtWidgets.QLabel("subsystem: ",self)
        subsys_comboBox_layout.addWidget(subsys_label)
        subsys_comboBox_layout.addWidget(self.subsys_comboBox)
        
        # 3rd from top
        simType_hLayout = QtWidgets.QHBoxLayout()        
        basis_label = QtWidgets.QLabel("basis: ",self)
        
        simType_hLayout.addWidget(basis_label)
        simType_hLayout.addWidget(self.basis_comboBox)
        
        # 4th from top
        cooper_N_hLayout = QtWidgets.QHBoxLayout()
        
        cooper_N_label = QtWidgets.QLabel(self)
        cooper_N_label.setText("cooper pairs max")
        
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
    
    def subsys_comboBox_changed_handler(self):
        self.subsystem_choice = self.subsys_comboBox.currentText()
        
    def basis_comboBox_changed_handler(self):
        self.simulation_basis = self.basis_comboBox.currentText()
    
    
    def fill_SL_names(self):
        self.SL_attributes_names = ["cooperN",
                                   "simulation_basis",
                                   "subsystem_choice"]
        
    
    def transfer_internal_to_widget(self):
        self.basis_comboBox.setCurrentText( self.simulation_basis )
        self.subsys_comboBox.setCurrentText( self.subsystem_choice ) 
        self.cooper_N_LineEdit.setText(str(self.cooperN))
            
        
        