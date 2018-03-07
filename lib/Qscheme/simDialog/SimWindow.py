from PyQt5 import QtCore,QtWidgets
from PyQt5.QtCore import Qt

from .SubscriptsWidget import SubscriptsWidget
from .SimulationRegimeWidget import SimulationRegimeWidget
from .ParametersSetupWidget import ParametersSetupWidget

import numpy as np
from collections import OrderedDict
from Qscheme.variables import Var

import pickle

from .SLBase import SLBase

class SimWindow(QtWidgets.QMainWindow,SLBase):
    def __init__(self,simulator):        
        super(SimWindow,self).__init__()
        
        self.setGeometry(100,100,500,300)
        self.setWindowTitle("Qcheme")
        
        self.simulator = simulator
        self.subscripts_widget = SubscriptsWidget(self)
        self.simulation_regime_widget = SimulationRegimeWidget(self)
        self.parameter_setup_widget = ParametersSetupWidget(self)        

        self.fill_SL_names()
        
        extractAction = QtWidgets.QAction("&Quit",self)
        extractAction.setShortcut("Ctrl+Q")
        extractAction.setStatusTip("Leave application")
        extractAction.triggered.connect(self.close)
        
        loadVars_file = QtWidgets.QAction("&open file",self)
        loadVars_file.setShortcut("Ctrl+O")
        loadVars_file.setStatusTip("file open dialog pops up")
        loadVars_file.triggered.connect(self.file_load)
        
        saveVars_file = QtWidgets.QAction("&save file",self)
        saveVars_file.setShortcut("Ctrl+S")
        saveVars_file.setStatusTip("file save dialog pops up")
        saveVars_file.triggered.connect(self.file_save)
        
        self.statusBar()
        
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu("&File")
        fileMenu.addAction(extractAction)
        
        saveLoadMenu = mainMenu.addMenu("&Save/Load")
        saveLoadMenu.addAction( loadVars_file )
        saveLoadMenu.addAction( saveVars_file )
        
        self.init_GUI()
        
    def init_GUI(self):
        '''Drawing GUI in the main frame'''
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.setAlignment(Qt.AlignTop)
        
        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        
        central_v_layout = QtWidgets.QVBoxLayout()
        central_widget.setLayout(v_layout)        
        
        central_line1_h_layout = QtWidgets.QHBoxLayout()
        
        central_line1_h_layout.addWidget(self.subscripts_widget)
        central_line1_h_layout.addWidget(self.simulation_regime_widget)
        central_v_layout.addLayout(central_line1_h_layout)        
        
        central_v_layout.addWidget(self.parameter_setup_widget)
        
        v_layout.addLayout(central_v_layout)
        
        bottom_buttons_layout = QtWidgets.QVBoxLayout()
        
        start_button = QtWidgets.QPushButton("Start simulation",self)
        start_button.clicked.connect(self.start_button_clicked_handler)
        bottom_buttons_layout.addWidget(start_button)
        
        v_layout.addLayout(bottom_buttons_layout)
        v_layout.addStretch(0)
        self.show()

    def start_button_clicked_handler(self):
        sim = self.simulator # name alias
        scheme = sim.scheme
        
        sim.cooperN = self.simulation_regime_widget.cooperN
        
        sweep_vars_intervals = OrderedDict()
        fixed_vars = []
        
        param_intervals = self.parameter_setup_widget.param_intervals

        mesh = OrderedDict()
        for param_sym in param_intervals:
            param_interval = param_intervals[param_sym]
            if( param_interval[0] != param_interval[1] ):
                sweep_vars_intervals[param_sym] = param_interval
                mesh[param_sym] = np.linspace(*param_interval)
            else:
                fixed_var = Var(param_sym,param_interval[0])
                fixed_vars.append( fixed_var )
                mesh[fixed_var.sym] = fixed_var.val
        
        #print("\nmesh: ",mesh)
        self.simulator.find_eigensystem_internal_params_product(mesh,4)
        
        spectr_idxs = [[0,1],[1,2],[2,3]]
        sweep_var = scheme.params[list(sweep_vars_intervals.keys())[0]]
        sim.plot2D_evals_from_var(sweep_var,
                                  fixed_vars,
                                  spectr_idxs)
     
    def fill_SL_names(self):
        self.SL_attributes_names = ["simulator"]
        self.SL_children_names = ["subscripts_widget",  
                                 "simulation_regime_widget", 
                                 "parameter_setup_widget"]
        
    def file_load(self):
        name = QtWidgets.QFileDialog.getOpenFileName(self,"Open File Name","Qscheme files","*.qsch")[0]
        with open(name,"rb") as file:
            load_dict = pickle.load(file)
        print(load_dict)
        self.load_from_dict_tree(load_dict)
        self.transfer_internal_to_widget_tree()
            
    def file_save(self):
        name = QtWidgets.QFileDialog.getSaveFileName(self,"Open File Name","","*.qsch")[0]          
        
        dump_dict = self.return_save_dict()
        with open(name,"wb") as file:
            pickle.dump(dump_dict,file)
    
    def printChildren(self,obj):
        #print( "children of ", obj )
        for child in obj.children():
            print( child )
            self.printChildren(child)
        #print( "end of ",obj )
        #print()
        
    def transfer_internal_to_widget(self):
        pass


