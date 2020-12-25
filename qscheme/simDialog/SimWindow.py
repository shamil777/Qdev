from PyQt5.QtWidgets import QMainWindow, QWidget
from PyQt5.QtWidgets import QVBoxLayout,QHBoxLayout
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QAction

from PyQt5.QtCore import Qt

from .SubscriptsWidget import SubscriptsWidget
from .SimulationRegimeWidget import SimulationRegimeWidget
from .ParametersSetupWidget import ParametersSetupWidget
from .SimulationDatasetsWidget import SimulationDatasetsWidget

from ..SLBase import SLBaseWidget
   

class SimWindow(QMainWindow,SLBaseWidget):        
    
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
        
        self.fill_SL_names()
        
        ## MAIN WINDOW WIDGETS SECTION START ##
        self.subscripts_widget = SubscriptsWidget(self)
        self.simulation_regime_widget = SimulationRegimeWidget(self)
        self.parameter_setup_widget = ParametersSetupWidget(self)
        self.simulation_datasets_widget = SimulationDatasetsWidget(self)
        ## MAIN WINDOW WIDGETS SECTION END
        
        # File menu item
        extractAction = QAction("&Quit",self)
        extractAction.setShortcut("Ctrl+Q")
        extractAction.setStatusTip("Leave application")
        extractAction.triggered.connect(self.close)
        
        # Savel/Load main menu items        
        loadVars_file = QAction("&open file",self)
        loadVars_file.setShortcut("Ctrl+O")
        loadVars_file.setStatusTip("file open dialog pops up")
        loadVars_file.triggered.connect(self.simulator.load_file_dialog)
        
        saveVars_file = QAction("&save file",self)
        saveVars_file.setShortcut("Ctrl+S")
        saveVars_file.setStatusTip("file save dialog pops up")
        saveVars_file.triggered.connect(self.simulator.save_file_dialog)
        
        # Import main menu item
        importPADs = QAction("import &netlist",self)
        importPADs.setStatusTip("import from .net netlist file format")
        importPADs.triggered.connect(self.simulator.import_netlist_handler)
        
        self.statusBar()
        
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu("&File")
        fileMenu.addAction(extractAction)
        
        saveLoadMenu = mainMenu.addMenu("&Save/Load")
        saveLoadMenu.addAction( loadVars_file )
        saveLoadMenu.addAction( saveVars_file )
        
        importMenu = mainMenu.addMenu("&Import")
        importMenu.addAction( importPADs )

        self.show()       
    
    def init_GUI(self):        
        '''Drawing GUI in the main frame'''  
        self.parameter_setup_widget.init_GUI()
        self.simulation_regime_widget.init_GUI()
        self.subscripts_widget.init_GUI()
        
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        
        central_v_layout = QVBoxLayout()
        self.central_widget.setLayout(central_v_layout)
        
        central_line1_h_layout = QHBoxLayout()
        
        central_line1_h_layout.addWidget(self.subscripts_widget)
        central_line1_h_layout.addWidget(self.simulation_regime_widget)
        central_v_layout.addLayout(central_line1_h_layout)        
        
        central_v_layout.addWidget(self.parameter_setup_widget)
        
        bottom_buttons_layout = QVBoxLayout()
        
        start_button = QPushButton("Start simulation",self)
        start_button.clicked.connect(self.start_button_clicked_handler)
        
        bottom_buttons_layout.addWidget(start_button)
        
        central_v_layout.addLayout(bottom_buttons_layout)
         
        self.gui_initialized = True

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
        self.SL_children_names = ["subscripts_widget",  
                                 "simulation_regime_widget", 
                                 "parameter_setup_widget"]
    
    def _printChildren(self,obj):
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


