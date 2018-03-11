import itertools

import numpy as np
import matplotlib.pyplot as plt

import networkx as nx
import pandas as pd

from collections import OrderedDict

from sympy.parsing.sympy_parser import parse_expr

from IPython.display import display

class VarSimKind:
    FIXED = "FIXED"
    SWEEP = "SWEEP"
    EQUATION = "EQUATION"

class SchemeSimulator:
    simulation_subsystems_keywords = ["Internal","External","Whole system"]
    simulation_basis_keywords = ["Node cooper pairs","Node phases"]
    
    def __init__(self, scheme, cooperN=None):
        self.scheme = scheme
        self.cooperN = cooperN
        
        simulation_datasets = pd.DataFrame(columns=["scheme_var_kinds","scheme_vars_settings",
                                                    "aux_var_kinds","aux_var_settings",
                                                    "simulation subsystem",
                                                    "simulation basis","simulation_basis_params"])
        
        
    
    def simulate(self,scheme_var_kinds,
                      scheme_var_settings,
                      aux_var_kinds,
                      aux_var_settings,
                      simulation_subsystem="Internal",
                      simulation_basis="Node cooper pairs",
                      **simulation_basis_params):
        errors = self.check_vars_sufficiency(scheme_var_kinds,scheme_var_settings,
                                    aux_var_kinds,aux_var_settings)
        print("errors:\n",errors)
        
    def check_vars_sufficiency(self,scheme_var_kinds,
                                    scheme_var_settings,
                                    aux_var_kinds,
                                    aux_var_settings):
        errors = OrderedDict() # {errorCode : reason of raising an error, ...}
        # Stage1
        for scheme_var_sym,var_setting in scheme_var_settings.items():
            if( scheme_var_kinds[scheme_var_sym] == VarSimKind.FIXED and 
                scheme_var_settings is None ):
                errors[scheme_var_sym] = "value is not set for fixed variable"
            elif( scheme_var_kinds[scheme_var_sym] == VarSimKind.SWEEP ):
                for setting in scheme_var_settings:
                    if( setting is None ):
                        errors[scheme_var_sym] = "value is not set for sweep variable"
            elif( scheme_var_kinds[scheme_var_sym] == VarSimKind.EQUATION ):
                if( scheme_var_settings is None ):
                    errors[scheme_var_sym] = "No equation is set for equation variable"
                    
        if( len(list(errors.keys())) > 0 ):
            return errors
        
        # Stage 2, parsing and checking
        graph_of_var_subs = nx.Graph()
        
        var_kinds = OrderedDict(**scheme_var_kinds,**aux_var_kinds)
        var_settings = OrderedDict(**scheme_var_settings,**aux_var_settings)
        parser_dict = OrderedDict([str(var_sym),var_sym] for var_sym in var_kinds.keys())
        for var_sym,var_kind in var_kinds.items():
            if(var_kind == VarSimKind.EQUATION):
                pars_result = parse_expr(var_settings[var_sym],local_dict=parser_dict)
                print(pars_result)
                
        return errors
    
    def find_eigensystem_internal(self,eigvals_N):
        self.scheme._construct_Hc_num_cooperN(self.cooperN)
        self.scheme._construct_Hj_num_cooperN(self.cooperN)
        H =  self.scheme.Hc_num_cooperN + self.scheme.Hj_num_cooperN
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
            
        return np.array(evals), einvects
    
    
    def find_eigensystem_on_mesh(self,params_list,mesh,eigvals_N):
        for point in mesh:
            print("finding eigval in point: '\n'{0} = {1}".format(params_list,point))
            self.scheme.assign_values_to_parameters(params_list,point)                
            #print( [(var.sym,var.val) for var in self.scheme.params.values()] )
            evals,einvects = self.find_eigensystem_internal(eigvals_N)
            
            self.points.append(self.scheme.get_params_values())
            self.einvects_list.append(einvects)
            self.evals_list.append(evals)
    
    def find_eigensystem_internal_product(self,params,eigvals_N):
        '''
        @description:
            Calculates eigensystem at every point in mesh
            that is constructed by outer product of params values lists
        @parameters:
            params - OrderedDict {"var1_sym":[var1_val1,var1_val2, ... ,var1_valN1],...}
        
        '''
        print(params)      
        error_symbols = self.check_vars_sufficiency(params)
        if( len(error_symbols) != 0 ):
            print("values for following symbols are not specified:" )
            print(error_symbols)
            return
        
        params = self.convert_single_to_lists(params)
        #print(params)
        
        vals_lists = [val for val in params.values()]
        params_list = list(params.keys())
        
        mesh = itertools.product(*vals_lists)

        self.find_eigensystem_on_mesh(params_list,mesh,eigvals_N)
            
        #print(self.evals_list[:,1]-self.evals_list[:,0])
    
    def find_eigensystem_internal_parametric(self,params,eigvals_N):
        '''
        @description:
            Calculates eigensystem at every point in mesh
            that is constructed params values lists
        @parameters:
            params - OrderedDict {"var1_sym":[var1_val1, ... ,var1_valN1],
                                  "var2_sym":[var1_val1, ... ,var1_valN1]}
                     Length of all lists should be equal, if not, lists which
                     length is lesser than maximal will be extended by repeating
                     its elemenets starting from the beggining.        
        '''
        error_symbols = self.check_vars_sufficiency(params)
        if( len(error_symbols) != 0 ):
            print("values for following symbols are not specified:" )
            print(error_symbols)
            return
        
        params = self.convert_single_to_lists(params)
        max_n_vals = max( list(map( len, params.values() )) )
        # filling vals_lists so, that all parameters lists
        # are having the same length.
        # With additional elements,
        # that are obtained by periodically repeating list elements
        # from the beggining of the list.
        vals_lists = []
        for sym,val in params.items():
            n_vals = len(val)
            if( n_vals < max_n_vals ):
                vals_lists.append([params[sym][i%n_vals] for i in range(max_n_vals)])
            elif( n_vals == max_n_vals ):
                vals_lists.append(params[sym])
            
        params_list = list(params.keys())
        mesh = zip(*vals_lists)
        print(mesh)
        self.find_eigensystem_on_mesh(params_list,mesh,eigvals_N)
        
    
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