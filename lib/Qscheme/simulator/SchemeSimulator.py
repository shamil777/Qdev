import itertools

import numpy as np
import matplotlib.pyplot as plt

import networkx as nx
import pandas as pd

from collections import OrderedDict

from .progressTimer import ProgressTimer

from sympy.parsing.sympy_parser import parse_expr

class VAR_KIND:
    FIXED = "FIXED"
    SWEEP = "SWEEP"
    EQUATION = "EQUATION"
    
class SIM_SUBSYS:
    INTERNAL = "Internal"
    COUPLING = "Coupling"
    WHOLE = "Whole"
    
class SIM_BASIS:
    COOPER = "Cooper"
    PHASE = "Phase"



class SchemeSimulator:
    simulation_subsystems_keywords = ["Internal","External","Whole system"]
    simulation_basis_keywords = ["Node cooper pairs","Node phases"]
    
    def __init__(self, scheme, cooperN=None, progress_window=None):
        self.progress_timer = ProgressTimer()
        
        self.scheme = scheme
        self.progress_window = progress_window
        self.cooperN = cooperN
        self.completed = 0
        
        simulation_datasets = pd.DataFrame(columns=["scheme_var_kinds","scheme_vars_settings",
                                                    "aux_var_kinds","aux_var_settings",
                                                    "simulation subsystem",
                                                    "simulation basis","simulation_basis_params","result"])
        self.parameters_consistent_graph = None
        self.leafs = None
        
        
    
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
        
        # result is stored as DataFrame object
        result = pd.DataFrame(columns=["E, GHz","eigen_vec"])
        
        # setting up union of the scheme and auxillary param dicts 
        # for data access convinience
        var_kinds = OrderedDict(**scheme_var_kinds,**aux_var_kinds)
        var_settings = OrderedDict(**scheme_var_settings,**aux_var_settings)

        # calculating the amount of points in parameters mesh
        iterations_n = 1
        for leaf_node_sym in self.leafs:
            if( var_kinds[leaf_node_sym] == VAR_KIND.SWEEP ):
                # multiplying by the number of points in sweep interval
                iterations_n *= var_settings[leaf_node_sym][2]
        
        # creating iterator that will return dict of the parameters with their new
        # values on each step of the mesh
        # starting from the leaf nodes
        
        # constructing leaf mesh
        leaf_var_settings = OrderedDict()
        for leaf_node_sym in self.leafs:
            if( var_kinds[leaf_node_sym] == VAR_KIND.SWEEP ):
                leaf_var_settings[leaf_node_sym] = np.linspace(*var_settings[leaf_node_sym])
            elif( var_kinds[leaf_node_sym] == VAR_KIND.FIXED ):
                leaf_var_settings[leaf_node_sym] = [var_settings[leaf_node_sym]]
        
        if( self.progress_window is not None ):
            self.progress_window.show()
            self.progress_window.update_progress()

        # local variables for time management
        self.dt_list = []

        # Cycling over leaf mesh and obtaining result.
        # Construction of dependent vars is made on the fly        
        leaf_mesh = itertools.product(*[itertools.product([key],val) for key,val in leaf_var_settings.items()])
        
        self.progress_timer.start(iterations_n)
        for iter_idx,leaf_point in enumerate(leaf_mesh):
            # getting the rest parameter values by evaluating
            # variables in equation tree
            vars_point = self.get_point_from_leaf_point(leaf_point,var_settings)
            
            # putting variables into scheme class
            for sym,var in self.scheme.params.items():
                var.val = vars_point[sym]
            
            # matrix diagonalization
            if( simulation_subsystem == SIM_SUBSYS.INTERNAL ):
                evals,evects = self.find_eigensystem_internal(eigenvals_num )
            elif( simulation_subsystem == SIM_SUBSYS.COUPLING ):
                raise NotImplementedError
            elif( simulation_subsystem == SIM_SUBSYS.WHOLE ):
                raise NotImplementedError
            
            # storing result
            for i in range(eigenvals_num):
                result.loc[i] = [evals[i],evects[i]]
            
            # updating time structure
            self.progress_timer.tick_dt()
                        
            # updating progress bar if necessary            
            if( self.progress_window is not None ):
                self.progress_window.update_progress(self.progress_timer)
            
            
            
        print(self.dt_list)         
        print(result)
        
        
        ## CALCULATION SECTION END ##
                
    def get_point_from_leaf_point(self,leaf_mesh_point,var_settings):
        vars_point = OrderedDict([(var.sym,var.val) for var in self.scheme.params.values()] )
        
        last_filled_nodes = self.leafs
        last_subs = OrderedDict([(sym,val) for sym,val in leaf_mesh_point])
        
        equation_graph = self.parameters_consistent_graph.copy()
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
                    
        self.parameters_consistent_graph = equations_dependency_graph
        self.leafs = leaf_nodes
        return errors
    
    def find_eigensystem_internal(self,eigvals_N):
        self.scheme._construct_Hc_num_cooperN(self.cooperN)
        self.scheme._construct_Hj_num_cooperN(self.cooperN)
        H =  self.scheme.Hc_num_cooperN + self.scheme.Hj_num_cooperN
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
            
        return np.array(evals), einvects
    
    
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

        
### USEFULL OLD CODE THAT CAN BE USED IN THE FUTURE ###

def array_eins_from_params_product(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,eigvals_N_pts):
    arr_evals = []
    arr_einvects = []
    
    fparams = [fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,eigvals_N_pts]
    for i,param in enumerate(fparams):
        if( not isinstance(param,list) and not isinstance(param,np.ndarray)):
            fparams[i] = [param]
            
    # sorted order is preserved
    prod = itertools.product(*fparams)
    
    for params in prod:
        f = params[0]*2*np.pi
        E_C = params[1]
        E_J = params[2]
        Csh = params[3]
        alpha = params[4]
        cooper_N = params[5]
        eigvals_N = params[6]
        
        H = H_J(E_J,alpha,f,cooper_N) + H_c(E_C,Csh,alpha,cooper_N)
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
        arr_evals.append(evals)
        arr_einvects.append(einvects)
        
    return arr_einvects, np.array(arr_evals)

def array_eins_from_params_lists(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,eigvals_N_pts):
    arr_evals = []
    arr_einvects = []
    # sorted order is preserved
    
    params = [fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,eigvals_N_pts]
    lengths = []
    for i, param in enumerate(params):
        if( isinstance(param,list) ):
            lengths.append(len(param))
        elif( isinstance(param,np.ndarray) ):
            lengths.append(param.shape[0])
        else:
            lengths.append(1)
            
    max_len = max(lengths)
    
    for i,param in enumerate(params):
        if( not isinstance(param,(list, np.ndarray)) ):
            params[i] = itertools.repeat(param,max_len)
        
    iter_list = zip( *params )
    
    for i,vals in enumerate(iter_list):
        f = vals[0]*2*np.pi
        E_C = vals[1]
        E_J = vals[2]
        Csh = vals[3]
        alpha = vals[4]
        cooper_N = vals[5]
        eigvals_N = vals[6]
        
        print('\r{:g} {:g} {:g} {:g} {:g} {:g} {:g}'.format(f,E_C,E_J,Csh,alpha,cooper_N,eigvals_N), end='')
        H = H_J(E_J,alpha,f,cooper_N) + H_c(E_C,Csh,alpha,cooper_N)
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
        arr_evals.append(evals)
        arr_einvects.append(einvects)
        
    return arr_einvects, np.array(arr_evals)

def array_coupling_from_parameters_product(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,multiplier):
    multiplier = 2*e*V_r0*beta/h # GHz
    
    g = []
    
    fparams = [fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,2]
    for i,param in enumerate(fparams):
        if( not isinstance(param,list) and not isinstance(param,np.ndarray)):
            fparams[i] = [param]
            
    # sorted order is preserved
    prod = itertools.product(*fparams)
    
    Evs, Engs = array_eins_from_params_product(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,2)
    Evs = np.array(Evs)
    
    for i,params in enumerate(prod):
        N = params[5]
        n1_op = qp.tensor(qp.charge(N),qp.identity(2*N+1))
        n2_op = qp.tensor(qp.identity(2*N+1),qp.charge(N))
        
        vec0 = Evs[i,0]
        vec1 = Evs[i,1]
        g.append( multiplier*(n1_op-n2_op).matrix_element(vec0.dag(),vec1) )
                 
    return np.array(g)

def array_coupling_from_parameters_lists(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,multiplier):
    multiplier = 2*e*V_r0*beta/h # GHz
    
    g = []
    
    Evs, Engs = array_eins_from_params_lists(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,2)
    Evs = np.array(Evs)
    
    params = [fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts,2]
    lengths = []
    for i, param in enumerate(params):
        if( isinstance(param,list) ):
            lengths.append(len(param))
        elif( isinstance(param,np.ndarray) ):
            lengths.append(param.shape[0])
        else:
            lengths.append(1)
            
    max_len = max(lengths)
    
    for i,param in enumerate(params):
        if( not isinstance(param,(list, np.ndarray)) ):
            params[i] = itertools.repeat(param,max_len)
        
    iter_list = zip( *params )
    
    for i, params in enumerate(iter_list):
        N = params[5]
        n1_op = qp.tensor(qp.charge(N),qp.identity(2*N+1))
        n2_op = qp.tensor(qp.identity(2*N+1),qp.charge(N))
        
        vec0 = Evs[i,0]
        vec1 = Evs[i,1]
        g.append( multiplier*(n1_op-n2_op).matrix_element(vec0.dag(),vec1) )
                 
    return np.array(g)

def array_eins_full_from_parameters_product(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts):
    pass