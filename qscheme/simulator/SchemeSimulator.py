import os
import sys
import csv
import pickle
import appdirs

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

from itertools import product
from collections import OrderedDict
from sympy.parsing.sympy_parser import parse_expr

from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QFileDialog

from .progressTimer import ProgressTimer
from ..SLBase import SLBase
from ..variables import Var
from ..Schematic.Scheme import Scheme
from ..simDialog import SimWindow
from ..simDialog.subWindows.SimulationOnGoing import SimulationOnGoingWindow
from .._KEYHASHABLE import SETTINGS_FIELD,SIM_SUBSYS,SIM_BASIS,VAR_KIND,IDX


class SchemeSimulator(SLBase):
    appData_dir = appdirs.user_data_dir()
    Qdev_folderName  = "Qdev"
    Qdev_settings_fileName = "Qdev_data.txt"
    Qdev_folder_path = os.path.join( appData_dir, Qdev_folderName )
    Qdev_settings_path = os.path.join( Qdev_folder_path, Qdev_settings_fileName )
    
    simulation_subsystems_keywords = ["Internal","External","Whole system"]
    simulation_basis_keywords = ["Node cooper pairs","Node phases"]
    
    def __init__(self, visualize=True, scheme=None ):
        SLBase.__init__(self)
        self.fill_SL_names()
        
        self.visualize = visualize
        
        self.progress_timer = None
        self.calculation_cancelled = False
        
        self.scheme = scheme
        self.cooperN = 1
        self.completed = 0
        
        self.simulation_datasets = pd.DataFrame(columns=["name","scheme_var_kinds","scheme_vars_settings",
                                                    "aux_var_kinds","aux_var_settings",
                                                    "simulation subsystem",
                                                    "simulation basis","simulation_basis_params","result"])

        self.result_last = None
        self.parameters_consistent_graph_current = None
        self.leafs_current = None
        
        ### GUI SECTION START ###
        if( not visualize ):
            pass
        else:
            if not QApplication.instance():
                self.app = QApplication(sys.argv)
            else:
                self.app = QApplication.instance()     

            self.main_window = SimWindow(self)
            self.progress_window = None
            self._reinit_progress_window()    
            
        # program settings, containing:
        # - info about last opened file
        # for more, see SettingsField static class
        self.settings = OrderedDict()
        self._load_settings() # loading settings if they are exist
        
        # if there were last session, then we load the appropriate file
        if( SETTINGS_FIELD.LAST_QDEV_FILE_OPENED in self.settings 
           and self.settings[SETTINGS_FIELD.LAST_QDEV_FILE_OPENED] is not None ):
            self.load_file(self.settings[SETTINGS_FIELD.LAST_QDEV_FILE_OPENED])            
        ### GUI SECTION END ###
        
    def fill_SL_names(self):
        self.SL_attributes_names = ["visualize",
                                    "scheme",
                                    "cooperN",
                                    "simulation_datasets"]
        self.SL_children_names = ["main_window"]
    
    def cancel_calculation(self):
        self.calculation_cancelled = True
    
    def find_eigenSystem(self,scheme_var_kinds,
                      scheme_var_settings,
                      aux_var_kinds,
                      aux_var_settings,
                      simulation_subsystem=SIM_SUBSYS.INTERNAL,
                      simulation_basis=SIM_BASIS.COOPER,
                      eigenvals_num = 5,
                      **simulation_basis_params):       
        
        # checking that parameters specifications
        # are self-sufficient
        errors = self.check_vars_sufficiency(scheme_var_kinds,scheme_var_settings,
                                    aux_var_kinds,aux_var_settings)
        print("errors:\n",errors)
        
        ## CALCULATION SECTION START ##
        
        # setting up some local variables
        self.cooperN = simulation_basis_params["cooper_N"]
        
        # setting up union of the scheme and auxillary param dicts 
        # for data access convinience
        var_kinds = OrderedDict(**scheme_var_kinds,**aux_var_kinds)
        var_settings = OrderedDict(**scheme_var_settings,**aux_var_settings)

        # calculating the amount of points in parameters mesh
        iterations_n = 1
        for leaf_node_sym in self.leafs_current:
            if( var_kinds[leaf_node_sym] == VAR_KIND.SWEEP ):
                # multiplying by the number of points in sweep interval
                iterations_n *= var_settings[leaf_node_sym][2]
        
        # creating iterator that will return dict of the parameters with their new
        # values on each step of the mesh
        # starting from the leaf nodes
        
        # constructing leaf mesh
        leaf_var_settings = OrderedDict()
        for leaf_node_sym in self.leafs_current:
            if( var_kinds[leaf_node_sym] == VAR_KIND.SWEEP ):
                leaf_var_settings[leaf_node_sym] = np.linspace(*var_settings[leaf_node_sym])
            elif( var_kinds[leaf_node_sym] == VAR_KIND.FIXED ):
                leaf_var_settings[leaf_node_sym] = [var_settings[leaf_node_sym]]
        
        # in case where a progress is present
        if( self.visualize is not None ):
            self.progress_window.show()
            self.progress_window.update_progress()

        # Cycling over leaf mesh and obtaining result.
        # Construction of dependent vars is made on the fly        
        leaf_mesh = product(*[product([key],val) for key,val in leaf_var_settings.items()])

        # making sure that calculation is not cancelled yet
        self.calculation_cancelled = False        
        
        # setting up timer for efficiency
        self.progress_timer = ProgressTimer()
        self.progress_timer.start(iterations_n)
        
        result = []
        if( self.main_window.parameter_setup_widget.sim_name is None ):
            pass # TODO: POP up window that shows that you ought to enter simulation name
            
        for iter_idx,leaf_point in enumerate(leaf_mesh):
            # getting the rest parameter values by evaluating
            # variables in equation tree
            vars_point = self.get_point_from_leaf_point(leaf_point,var_settings)
            
            # putting variables into scheme class
            for sym,var in self.scheme.params.items():
                var.val = vars_point[sym]
            
            # matrix diagonalization
            if( simulation_subsystem == SIM_SUBSYS.INTERNAL ):
                evals, evects = self.find_eigensystem_internal( eigenvals_num )
            elif( simulation_subsystem == SIM_SUBSYS.COUPLING ):
                raise NotImplementedError
            elif( simulation_subsystem == SIM_SUBSYS.WHOLE ):
                raise NotImplementedError
                
            # storing spectra of the current mesh point
            result.append( [vars_point,evals,evects] ) # storing point result to list
            
            # updating time structure
            self.progress_timer.tick_dt()
                        
            # updating progress bar if necessary            
            if( self.progress_window is not None ):
                self.progress_window.update_progress(self.progress_timer)
            
            # whether the calculation is cancelled
            if( self.calculation_cancelled is True ):
                break
            
        self.progress_window.plot_btn_activate()
        result = np.array(result)
        name = self.main_window.parameter_setup_widget.sim_name
        res_dict = {"scheme_var_kinds":scheme_var_kinds,
                    "simulation_name":name,
                    "scheme_vars_settings":scheme_var_settings,
                    "aux_var_kinds":aux_var_kinds,
                    "aux_var_settings":aux_var_settings,
                    "simulation subsystem":simulation_subsystem,
                    "simulation basis":simulation_basis,
                    "simulation_basis_params":simulation_basis_params,
                    "result":result}
        self.simulation_datasets = self.simulation_datasets.append(res_dict,ignore_index=True)  
        
        ## CALCULATION SECTION END ##
    
    def _print_array_info( self, name, array, data=False ):
        print(name + ": ",type(array),"| data type: ",array.dtype)
        print("shape: ", array.shape)
        if( data ):
            print(array,'\n')
            print(array.T)
    
    def plot_last_result(self):
        last_res = self.simulation_datasets["result"].iloc[-1:].values[0]
        sweep_var = Var("t")
        x = np.array([pt[sweep_var.sym] for pt in last_res[:,IDX.POINT]])        
        y_list = last_res[:,IDX.EVALS]
        
        y_list = np.vstack(y_list).T
        
        spectr_idxs = [[0,1],[1,2],[2,3]]
        
        for i,idx in enumerate(spectr_idxs):
            plt.plot(x,y_list[idx[1]] - y_list[idx[0]], label="$E_" + str(idx[1]) + " - E_"+ str(idx[0]) + "$")
            
        # plot cosmetics
        plt.ylabel(r"E, GHz")
        plt.xlabel(r"$" + str(sweep_var.sym) + "$")
        plt.legend()
        plt.grid()
        
        plt.show()
            
            
    
    def get_point_from_leaf_point(self,leaf_mesh_point,var_settings):
        vars_point = OrderedDict([(var.sym,var.val) for var in self.scheme.params.values()] )
        
        last_filled_nodes = self.leafs_current
        last_subs = OrderedDict([(sym,val) for sym,val in leaf_mesh_point])
        
        equation_graph = self.parameters_consistent_graph_current.copy()
        while( len(last_filled_nodes) > 0 ):
            new_filled_nodes = []
            new_subs = OrderedDict()
            
            vars_point.update(last_subs)
            
            # substitution of last_filled_nodes into their ancestors
            for in_edge in equation_graph.in_edges(last_filled_nodes):
                node_sym = in_edge[0]
                # updating equation in the ancestor node
                parsed_eq = equation_graph.nodes[node_sym]['parsed_eq']
                equation_graph.nodes[node_sym]['parsed_eq'] = parsed_eq.subs(last_subs)
                    
                # if this node is fully numeric we add it as a source to the
                # next algorithm step
                try:
                    new_val = float(equation_graph.nodes[node_sym]['parsed_eq'])
                except ValueError:
                    continue
                else:
                    new_filled_nodes.append(node_sym)
                    new_subs[node_sym] = new_val
            
            last_filled_nodes = new_filled_nodes
            last_subs = new_subs
            
        return vars_point
                
            
        
        
    def check_vars_sufficiency(self,scheme_var_kinds,
                                    scheme_var_settings,
                                    aux_var_kinds,
                                    aux_var_settings):
        errors = OrderedDict() # {errorCode : reason of raising an error, ...}
        # Stage1
        for scheme_var_sym,var_setting in scheme_var_settings.items():
            if( scheme_var_kinds[scheme_var_sym] == VAR_KIND.FIXED and 
                scheme_var_settings is None ):
                errors[scheme_var_sym] = "value is not set for fixed variable"
            elif( scheme_var_kinds[scheme_var_sym] == VAR_KIND.SWEEP ):
                for setting in scheme_var_settings:
                    if( setting is None ):
                        errors[scheme_var_sym] = "value is not set for sweep variable"
            elif( scheme_var_kinds[scheme_var_sym] == VAR_KIND.EQUATION ):
                if( scheme_var_settings is None ):
                    errors[scheme_var_sym] = "No equation is set for equation variable"
                    
        if( len(list(errors.keys())) > 0 ):
            return errors
        
        # Stage 2, parsing and checking
        equations_dependency_graph = nx.MultiDiGraph()
        
        var_kinds = OrderedDict(**scheme_var_kinds,**aux_var_kinds)
        var_settings = OrderedDict(**scheme_var_settings,**aux_var_settings)
        parser_dict = OrderedDict([str(var_sym),var_sym] for var_sym in var_kinds.keys())
        
        # constructing equation depencies graph
        for var_sym,var_kind_str in var_kinds.items():
            if(var_kind_str == VAR_KIND.EQUATION):
                pars_result = parse_expr(var_settings[var_sym],local_dict=parser_dict)
                equations_dependency_graph.add_node(var_sym,var_kind=var_kind_str,parsed_eq=pars_result)
                for eq_sym in pars_result.free_symbols:
                    equations_dependency_graph.add_node(eq_sym,var_kind=var_kinds[eq_sym])
                    equations_dependency_graph.add_edge(var_sym,eq_sym)
            elif( var_kind_str == VAR_KIND.FIXED or var_kind_str == VAR_KIND.SWEEP ):
                equations_dependency_graph.add_node(var_sym,var_kind=var_kind_str)
            
        
        # if graph contains any cycles, we terminate the simulation
        try:
            cycles = nx.find_cycle(equations_dependency_graph)
        except nx.NetworkXNoCycle:
            pass
        else:
            errors["DependencyError"] = "Cycle is found in equation depencies graph"
            return errors
        
        # if graph leafs is not of type 'FIXED' or 'SWEEP'
        # terminate simulation
        
        leaf_nodes = []
        
        for node_var_sym,node_var_kind in equations_dependency_graph.nodes.data('var_kind'):
            successors = equations_dependency_graph.successors(node_var_sym)
            # in case this node has no successors
            if( len(list(successors)) == 0 ):
                if( node_var_kind == VAR_KIND.EQUATION ):
                    errors["DependencyError"] = "Dependency graph contains no cycles, \
                                                 but one or more of it's leafs has type 'EQUATION'"
                    return errors
                else:
                    leaf_nodes.append(node_var_sym)    
                    
        self.parameters_consistent_graph_current = equations_dependency_graph
        self.leafs_current = leaf_nodes
        return errors
    
    def find_eigensystem_internal(self,eigvals_N):
        self.scheme._construct_Hc_num_cooperN(self.cooperN)
        self.scheme._construct_Hj_num_cooperN(self.cooperN)
        H =  self.scheme.Hc_num_cooperN + self.scheme.Hj_num_cooperN
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
        einvects = np.array([vec.data for vec in einvects])
        return evals, einvects
    
    
    def find_eigensystem_on_mesh(self,params_list,mesh,eigvals_N):
        for point in mesh:
            self.scheme.assign_values_to_parameters(params_list,point)                
            #print( [(var.sym,var.val) for var in self.scheme.params.values()] )
            evals,einvects = self.find_eigensystem_internal(eigvals_N)
            
            self.points.append(self.scheme.get_params_values())
            self.einvects_list.append(einvects)
            self.evals_list.append(evals)
        
    
    def convert_single_to_lists(self,params):
        for param_sym,param_val in params.items():
            if( not isinstance(param_val,list) and not isinstance(param_val,np.ndarray)):
                params[param_sym] = [param_val]
                
        return params
                
        
    def plot2D_evals_from_var( self, sweep_var, fixed_vars, spectr_idxs ):
        '''
        @description:
        @parameters:
            sweep_var - variable along the x axis of the plot
            fixed_vars - list with fixed variables
        '''
        params_keys = list(self.scheme.params.keys())
        
        # collect all points corresponding to fixed vars
        # Optimized in way that all array checking is performed
        # inside the numpy
        
        # constructing byte-mask for points
        fixed_vars_mask = np.ones(len(self.scheme.params), dtype=np.float64)
        
        sweep_var_idx = list(params_keys).index(sweep_var.sym)
        fixed_vars_mask[sweep_var_idx] = 0
        
        # rebuilding fixed_vars into point-compatible format
        fixed_vars_as_point = np.zeros(len(fixed_vars_mask), dtype=np.float64)
        
        # collecting points of interest, by iterating through
        # all data and checking
        for i,fixed_var in enumerate(fixed_vars):
            fixed_param_true_idx = params_keys.index(fixed_var.sym)
            fixed_vars_as_point[fixed_param_true_idx] += fixed_vars[i].val
            
        # collecting points of interest indexes
        points_idxs = []
        for i,point in enumerate(self.points):
            if( (point*fixed_vars_mask - fixed_vars_as_point).all() == 0 ):
                points_idxs.append(i)
        
        # gathering x,y data
        self.points = np.array(self.points)
        x = self.points[points_idxs,sweep_var_idx]
        self.points = list(self.points)
        
        # gathering y-data
        self.einvects_list = np.array(self.einvects_list)
        self.evals_list = np.array(self.evals_list)
        y_list = []*len(spectr_idxs)
        for idx in spectr_idxs:
            # y_list contains E_idx - E_0 for all idx in specr_idx
            y_list.append(self.evals_list[points_idxs,idx[1]] - self.evals_list[points_idxs,idx[0]])
        self.einvects_list = list(self.einvects_list)
        self.evals_list = list(self.evals_list)
        
        # plotting
        for idx,y in zip(spectr_idxs,y_list):
            plt.plot(x,y, label="$E_" + str(idx[1]) + " - E_"+ str(idx[0]) + "$")
            
        # plot cosmetics
        plt.ylabel(r"E, GHz")
        plt.xlabel(r"$" + str(sweep_var.sym) + "$")
        plt.legend()
        plt.grid()
        
        plt.show()
        
    ###### API SECTION START ######
    ## GUI SECTION START ##
    def GUI_loop(self):
        '''
        @description: launches application main window
        '''
        return self.app.exec_()
        
    def _reinit_progress_window(self):
        if( not self.visualize ):
            return
        
        if( self.progress_window is not None ):
            self.progress_window.setParent(None)
            
        self.progress_window = SimulationOnGoingWindow(plot_func_ref=self.plot_last_result)
    ## GUI SECTION END ##
    
    ## SCHEME LOAD SECTION START ##
    def import_netlist_handler(self):
        file_name = QFileDialog.getOpenFileName(self.main_window,"Open File Name","netlist files","*.net")[0]
        self.scheme = Scheme(file_path=file_name)
        self.main_window.scheme = self.scheme
        if( not self.main_window.gui_initialized ):
            self.main_window.init_GUI()
        else:
            self.main_window.parameter_setup_widget.scheme_subscripts_changed_handler()
    ## SCHEME LOAD SECTION END ##
    
    ## FILE SAVE SECTION START ##        
    def save_file(self,filepath):
        '''
        @description: saves this simulator class state to drive
        '''
        dump_dict = self.return_save_dict()
        self._preSave_unpicklable_refs_nullifier()
        try:
            with open(filepath,"wb") as file:
                pickle.dump(dump_dict,file)
        except FileNotFoundError:
            return
            
        self.settings[SETTINGS_FIELD.LAST_QDEV_FILE_OPENED] = filepath
        self._save_settings()
        self._unpickable_refs_restorator() # to continue working with all links restored
    
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
        pass
        
    def save_file_dialog(self):
        '''
        @description:   shows save file dialog and then invokes
                        save_state_to_file(filename) with filename
                        derived from dialog.
        @return: None
        '''
        filepath = QFileDialog.getSaveFileName(self.main_window,"Open File Name","","*.qsch")[0]          
        self.save_file(filepath)
    ## FILE SAVE SECTION END ##
        
    ## FILE LOAD SECTION START ##
    def load_file(self,filepath):
        with open(filepath,"rb") as file:
            load_dict = pickle.load(file)

        self.load_from_dict_tree(load_dict)
        self.handle_loadData_tree()
        self._unpickable_refs_restorator() 

        if( not self.main_window.gui_initialized ):            
            self.main_window.init_GUI()
        
        self.settings[SETTINGS_FIELD.LAST_QDEV_FILE_OPENED] = filepath
        self._save_settings()
        
    def _unpickable_refs_restorator(self):
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
        pass
        
    def load_file_dialog(self):
        filepath = QFileDialog.getOpenFileName(self.main_window,"Open File Name","qscheme files","*.qsch")[0]
        self.load_file(filepath)
    ## FILE LOAD SECTION START ##
    
    ## PROGRAM DATA ON DRIVE SECTION START ##
    def _load_settings(self):
        # if there is no such file or dir -> return
        if( not os.path.exists(self.Qdev_settings_path) ):
            return
        
        with open(self.Qdev_settings_path,"r") as file:
            reader = csv.reader(file,delimiter='=')
            for row in reader:
                if( len(row) > 1 ):
                    self.settings[row[0]] = row[1]
        
    def _save_settings(self):
        # if there is no such dir, create it
        if( not os.path.isdir(self.Qdev_folder_path) ):
            os.mkdir(self.Qdev_folder_path)
            
        # open file, overriding previous content
        with open(self.Qdev_settings_path,"w") as file:
            writer = csv.writer(file,delimiter='=')
            for key,val in self.settings.items():
                writer.writerow([key,val])
    ## PROGRAM DATA ON DRIVE SECTION END ##      
    ###### API SECTION END ######