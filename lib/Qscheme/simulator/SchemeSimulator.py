import itertools

import numpy as np
import matplotlib.pyplot as plt

class SchemeSimulator:
    def __init__(self, scheme, cooperN):
        self.scheme = scheme
        self.cooperN = cooperN
        
        self.points = []
        self.einvects_list =[]
        self.evals_list = []
    
    def find_eigensystem_internal(self,eigvals_N):
        self.scheme._construct_Hc_num_cooperN(self.cooperN)
        self.scheme._construct_Hj_num_cooperN(self.cooperN)
        H =  self.scheme.Hc_num_cooperN + self.scheme.Hj_num_cooperN
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
            
        return np.array(evals), einvects
    
    
    def find_eigensystem_internal_params_product(self,params,eigvals_N):
        '''
        @description:
            Calculates eigensystem at every point in mesh
            that is constructed by outer product of params values lists
        @parameters:
            params - OrderedDict {"var1_sym":[var1_val1,var1_val2, ... ,var1_valN1],...}
        
        '''
        
        
        error_symbols = self.check_vars_sufficiency(params)
        if( len(error_symbols) != 0 ):
            print("values for following symbols are not specified:" )
            print(error_symbols)
            return
        
        params = self.convert_single_to_lists(params)
        print(params)
        
        vals_lists = [val for val in params.values()]
        params_list = list(params.keys())
        mesh = itertools.product(*vals_lists)
        
        for point in mesh:
            self.scheme.assign_values_to_parameters(params_list,point)                
            print( [(var.sym,var.val) for var in self.scheme.params.values()] )
            evals,einvects = self.find_eigensystem_internal(eigvals_N)
            
            self.points.append(self.scheme.get_params_values())
            self.einvects_list.append(einvects)
            self.evals_list.append(evals)
            
        self.points = np.array(self.points)
        self.einvects_list = np.array(self.einvects_list)
        self.evals_list = np.array(self.evals_list)
        print(self.evals_list[:,1]-self.evals_list[:,0])
        
    def convert_single_to_lists(self,params):
        for param_sym,param_val in params.items():
            if( not isinstance(param_val,list) ):
                params[param_sym] = [param_val]
                
        return params
                
    def check_vars_sufficiency(self,params):
        syms_list = [sym for sym in params]
        error_syms = []
        for var in self.scheme.params.values():
            if( var.sym in syms_list ):
                continue
            elif( var.val is None ):
                error_syms.append(var.sym)
        
        return error_syms
    
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
        x = self.points[points_idxs,sweep_var_idx]
        
        # gathering y-data
        y_list = []*len(spectr_idxs)
        for idx in spectr_idxs:
            # y_list contains E_idx - E_0 for all idx in specr_idx
            y_list.append(self.evals_list[points_idxs,idx[1]] - self.evals_list[points_idxs,idx[0]])
        
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