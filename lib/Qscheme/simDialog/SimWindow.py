from PyQt5 import QtCore,QtWidgets
from PyQt5.QtCore import Qt

from .SubscriptsWidget import SubscriptsWidget
from .SimulationRegimeWidget import SimulationRegimeWidget
from .ParametersSetupWidget import ParametersSetupWidget
from .subWindows.SimulationOnGoing import SimulationOnGoingWindow

from collections import OrderedDict

from Qscheme.schematic.Scheme import Scheme
from Qscheme.simulator import SchemeSimulator

import pickle

from .SLBase import SLBase

import appdirs
import os
import csv

class SettingsField():
    LAST_QDEV_FILE_OPENED = "LAST_QDEV_FILE_OPENED"
   

class SimWindow(QtWidgets.QMainWindow,SLBase):
    appData_dir = appdirs.user_data_dir()
    Qdev_folderName  = "Qdev"
    Qdev_settings_fileName = "Qdev_data.txt"
    Qdev_folder_path = os.path.join( appData_dir, Qdev_folderName )
    Qdev_settings_path = os.path.join( Qdev_folder_path, Qdev_settings_fileName )        
    
    def __init__(self,simulator=None):        
        super(SimWindow,self).__init__( )

        self.gui_initialized = False
        
        self.setGeometry(100,100,500,300)
        self.setWindowTitle("Qcheme")

        self.simulator = simulator
        self.scheme = None
        if( self.simulator is not None and 
            self.simulator.scheme is not None ):
            self.scheme = self.simulator.scheme    

        ## MAIN WINDOW WIDGETS SECTION START ##
        self.subscripts_widget = SubscriptsWidget(self)
        self.simulation_regime_widget = SimulationRegimeWidget(self)
        self.parameter_setup_widget = ParametersSetupWidget(self) 
        ## MAIN WINDOW WIDGETS SECTION END

        ## ADDITIONAL WINDOWS SECTION START ##
        self.progress_window = SimulationOnGoingWindow(self)
        ## ADDITIONAL WINDOWS SECTION END ##  
        
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
        loadVars_file.triggered.connect(self.file_load_dialog)
        
        saveVars_file = QtWidgets.QAction("&save file",self)
        saveVars_file.setShortcut("Ctrl+S")
        saveVars_file.setStatusTip("file save dialog pops up")
        saveVars_file.triggered.connect(self.file_save_dialog)
        
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
        
        # program settings, containing:
        # - info about last opened file
        # for more, see SettingsField static class
        self.settings = OrderedDict()
        self._load_settings() # loading settings if they are exist
        
        # if there were last session, then we load the appropriate file
        if( SettingsField.LAST_QDEV_FILE_OPENED in self.settings 
           and self.settings[SettingsField.LAST_QDEV_FILE_OPENED] is not None ):
            self.file_load(self.settings[SettingsField.LAST_QDEV_FILE_OPENED])            
        elif( self.simulator is not None and self.scheme is not None ):
            self.init_GUI()
            
        self.show()
    
    def _load_settings(self):
        # if there is no such file or dir -> return
        if( not os.path.exists(SimWindow.Qdev_settings_path) ):
            return
        
        with open(SimWindow.Qdev_settings_path,"r") as file:
            reader = csv.reader(file,delimiter='=')
            for row in reader:
                if( len(row) > 1 ):
                    self.settings[row[0]] = row[1]
        
    def _save_settings(self):
        # if there is no such dir, create it
        if( not os.path.isdir(SimWindow.Qdev_folder_path) ):
            os.mkdir(SimWindow.Qdev_folder_path)
            
        # open file, overriding previous content
        with open(SimWindow.Qdev_settings_path,"w") as file:
            writer = csv.writer(file,delimiter='=')
            for key,val in self.settings.items():
                writer.writerow([key,val])
            
    
    def init_GUI(self):
        '''Drawing GUI in the main frame'''  
        self.parameter_setup_widget.init_GUI()
        self.simulation_regime_widget.init_GUI()
        self.subscripts_widget.init_GUI()
        
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
        visualize_button = QtWidgets.QPushButton("Visualize",self)
        visualize_button.clicked.connect(self.visualize_button_clicked_handler)
        
        bottom_buttons_layout.addWidget(start_button)
        bottom_buttons_layout.addWidget(visualize_button)
        
        v_layout.addLayout(bottom_buttons_layout)

        self.simulator.progress_window = self.progress_window
        
        self.gui_initialized = True
        #print("children: ", [child.__class__ for child in self.children()])

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
    
    def visualize_button_clicked_handler(self):
        self.simulator.plot_last_result()
    
    def fill_SL_names(self):
        self.SL_attributes_names = ["simulator"]
        self.SL_children_names = ["subscripts_widget",  
                                 "simulation_regime_widget", 
                                 "parameter_setup_widget"]
    
    def _preSave_unpicklable_refs_nullifier(self):
        '''
        @description: ## FOR INTERNAL USAGE ##
                Pre-save unserializable references nullifier.
                Pickle.dump will return error if some of the object are containing
            unserializable objects, like "simulator.progress_window" is widget and
            hence not serializable.
                Used in pair with "_afterLoad_unpickable_refs_restorator"
            that restores all nullified references if possible.
        @args: None
        @return: None
        '''
        self.simulator.progress_window = None
    
    def file_save_dialog(self):
        name = QtWidgets.QFileDialog.getSaveFileName(self,"Open File Name","","*.qsch")[0]          
        self.file_save(name)
        
    def file_save(self,name):
        dump_dict = self.return_save_dict()
        self._preSave_unpicklable_refs_nullifier()
        with open(name,"wb") as file:
            pickle.dump(dump_dict,file)
            
        self.settings[SettingsField.LAST_QDEV_FILE_OPENED] = name
        self._save_settings()
    
    def _afterLoad_unpickable_refs_restorator(self):
        '''
        @description: ## FOR INTERNAL USAGE ##
                After-load unserializable references restoration.
                Pickle.dump will return error if some of the object are containing
            unserializable objects, like "simulator.progress_window" is widget and
            hence not serializable. This function restores the unpicklable links
            after load process is over.
                Used in pair with "_preSave_unpicklable_refs_nullifier"
            that nullifies corresponding references (set them to None).
        @args: None
        @return: None
        '''
        self.simulator.progress_window = self.progress_window
        
    def file_load_dialog(self):
        name = QtWidgets.QFileDialog.getOpenFileName(self,"Open File Name","Qscheme files","*.qsch")[0]
        self.file_load(name)
        
    def file_load(self,name):
        with open(name,"rb") as file:
            load_dict = pickle.load(file)

        if( not self.gui_initialized ):
            self.simulator = load_dict["simulator"]
            self.scheme = self.simulator.scheme
        
            self.init_GUI()
            
        self.load_from_dict_tree(load_dict)
        self.transfer_internal_to_widget_tree()
        self._afterLoad_unpickable_refs_restorator() 
        
        self.settings[SettingsField.LAST_QDEV_FILE_OPENED] = name
        self._save_settings()
            
    
    def import_netlist_handler(self):
        file_name = QtWidgets.QFileDialog.getOpenFileName(self,"Open File Name","netlist files","*.net")[0]
        self.scheme = Scheme(file_path=file_name)
        self.simulator = SchemeSimulator(self.scheme)
        if( not self.gui_initialized ):
            self.init_GUI()
        else:
            self.parameter_setup_widget.scheme_subscripts_changed_handler()
    
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


