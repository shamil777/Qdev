from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from .SLBase import SLBase
from ..simulator import SchemeSimulator

class SimulationRegimeWidget(QtWidgets.QWidget,SLBase):
    
    def __init__(self, parent=None, flags=Qt.WindowFlags()):
        super(SimulationRegimeWidget,self).__init__(parent,flags)
        
        self.simulation_basis_strs = SchemeSimulator.simulation_basis_keywords
        self.subsystem_choice_strs = SchemeSimulator.simulation_subsystems_keywords
        
        # Save/Load variables section START #    
        self.cooperN = 5
        self.simulation_basis = None
        self.subsystem_choice = None
        # Save/Load variables section END #
        
        self.init_GUI()
        self.fill_SL_names()
        
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
        self.subsys_comboBox = QtWidgets.QComboBox()
        self.subsys_comboBox.addItems(self.subsystem_choice_strs)
        self.subsys_comboBox.currentIndexChanged.connect(self.subsys_comboBox_changed_handler)
        self.subsys_comboBox_changed_handler()
        
        subsys_comboBox_layout.addWidget(subsys_label)
        subsys_comboBox_layout.addWidget(self.subsys_comboBox)
        
        # 3rd from top
        simType_hLayout = QtWidgets.QHBoxLayout()
        
        basis_label = QtWidgets.QLabel("basis: ",self)
        self.basis_comboBox = QtWidgets.QComboBox(self)
        self.basis_comboBox.addItems(self.simulation_basis_strs)
        self.basis_comboBox.currentIndexChanged.connect(self.basis_comboBox_changed_handler)
        self.basis_comboBox_changed_handler()
        
        simType_hLayout.addWidget(basis_label)
        simType_hLayout.addWidget(self.basis_comboBox)
        
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
    
    def subsys_comboBox_changed_handler(self):
        self.subsystem_choice = self.subsys_comboBox.currentText()
        
    def basis_comboBox_changed_handler(self):
        self.simulation_basis = self.basis_comboBox.currentText()
    
    
    def fill_SL_names(self):
        self.SL_attributes_names = ["cooperN",
                                   "simulation_basis",
                                   "subsystem_choice"]
        
    
    def transfer_internal_to_widget(self):
        print(self.simulation_basis)
        print(self.subsystem_choice)
        self.basis_comboBox.setCurrentIndex(
                self.simulation_basis_strs.index(self.simulation_basis)
                )
        self.subsys_comboBox.setCurrentIndex(
                self.subsystem_choice_strs.index(self.subsystem_choice))
        self.cooper_N_LineEdit.setText(str(self.cooperN))
            
        
        