from PyQt5 import QtCore,QtWidgets
from PyQt5.QtCore import Qt

from .SubscriptsWidget import SubscriptsWidget
from .SimulationRegimeWidget import SimulationRegimeWidget
from .ParametersSetupWidget import ParametersSetupWidget
from .subWindows.SimulationOnGoing import SimulationOnGoingWindow


from Qscheme.schematic.Scheme import Scheme
from Qscheme.simulator import SchemeSimulator

import pickle

from .SLBase import SLBase

class SimWindow(QtWidgets.QMainWindow,SLBase):
    def __init__(self,simulator=None):        
        super(SimWindow,self).__init__()
        
        self.gui_initialized = False
        
        self.setGeometry(100,100,500,300)
        self.setWindowTitle("Qcheme")
        
        self.simulator = simulator
        self.scheme = None
        if( self.simulator is not None and 
            self.simulator.scheme is not None ):
            self.scheme = self.simulator.scheme    
        
        self.progress_window = SimulationOnGoingWindow(self)
 
        self.fill_SL_names()
        
        # File menu item
        extractAction = QtWidgets.QAction("&Quit",self)
        extractAction.setShortcut("Ctrl+Q")
        extractAction.setStatusTip("Leave application")
        extractAction.triggered.connect(self.close)
        
        # Savel/Load main menu items        
        loadVars_file = QtWidgets.QAction("&open file",self)
        loadVars_file.setShortcut("Ctrl+O")
        loadVars_file.setStatusTip("file open dialog pops up")
        loadVars_file.triggered.connect(self.file_load)
        
        saveVars_file = QtWidgets.QAction("&save file",self)
        saveVars_file.setShortcut("Ctrl+S")
        saveVars_file.setStatusTip("file save dialog pops up")
        saveVars_file.triggered.connect(self.file_save)
        
        # Import main menu item
        importPADs = QtWidgets.QAction("import &netlist",self)
        importPADs.setStatusTip("import from .net netlist file format")
        importPADs.triggered.connect(self.import_netlist_handler)
        
        self.statusBar()
        
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu("&File")
        fileMenu.addAction(extractAction)
        
        saveLoadMenu = mainMenu.addMenu("&Save/Load")
        saveLoadMenu.addAction( loadVars_file )
        saveLoadMenu.addAction( saveVars_file )
        
        importMenu = mainMenu.addMenu("&Import")
        importMenu.addAction( importPADs )
        
        if( self.simulator is not None and self.scheme is not None ):
            self.init_GUI()
            
        self.show()
        
    def init_GUI(self):
        '''Drawing GUI in the main frame'''
        self.subscripts_widget = SubscriptsWidget(self)
        self.simulation_regime_widget = SimulationRegimeWidget(self)
        self.parameter_setup_widget = ParametersSetupWidget(self) 
        
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.setAlignment(Qt.AlignTop)
        
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        
        central_v_layout = QtWidgets.QVBoxLayout()
        self.central_widget.setLayout(v_layout)        
        
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

        self.simulator.progress_window = self.progress_window
        
        self.gui_initialized = True
        print("children: ", [child.__class__ for child in self.children()])

    def start_button_clicked_handler(self):
        # name alias for shorter reference to objects
        scheme_var_grid = self.parameter_setup_widget.scheme_vars_grid
        aux_var_grid = self.parameter_setup_widget.aux_vars_grid
        
        self.simulator.find_eigenSystem(scheme_var_grid.var_kinds,
                                scheme_var_grid.var_settings,
                                aux_var_grid.var_kinds,
                                aux_var_grid.var_settings,
                                simulation_subsystem=self.simulation_regime_widget.subsystem_choice,
                                simulation_basis=self.simulation_regime_widget.simulation_basis,
                                cooper_N=self.simulation_regime_widget.cooperN)
     
    def fill_SL_names(self):
        self.SL_attributes_names = ["simulator"]
        self.SL_children_names = ["subscripts_widget",  
                                 "simulation_regime_widget", 
                                 "parameter_setup_widget"]
        
    def file_load(self):
        name = QtWidgets.QFileDialog.getOpenFileName(self,"Open File Name","Qscheme files","*.qsch")[0]
        with open(name,"rb") as file:
            load_dict = pickle.load(file)
        
        if( not self.gui_initialized ):
            self.simulator = load_dict["simulator"]
            self.scheme = self.simulator.scheme
        
            self.init_GUI()
        
        self.load_from_dict_tree(load_dict)
        self.transfer_internal_to_widget_tree()
            
    def file_save(self):
        name = QtWidgets.QFileDialog.getSaveFileName(self,"Open File Name","","*.qsch")[0]          
        
        dump_dict = self.return_save_dict()
        with open(name,"wb") as file:
            pickle.dump(dump_dict,file)
    
    def import_netlist_handler(self):
        file_name = QtWidgets.QFileDialog.getOpenFileName(self,"Open File Name","netlist files","*.net")[0]
        self.scheme = Scheme(file_path=file_name)
        self.simulator = SchemeSimulator(self.scheme)
        if( not self.gui_initialized ):
            self.init_GUI()
        else:
            self.parameter_setup_widget.scheme_subscripts_changed_handler()
    
    def printChildren(self,obj):
        #print( "children of ", obj )
        for child in obj.children():
            print( child )
            self.printChildren(child)
        #print( "end of ",obj )
        #print()
        
    def transfer_internal_to_widget(self):
        pass
    
    def fit_window_to_content(self):
        self.resize(self.minimumSizeHint())


