# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 20:24:02 2017

@author: user
"""
#from sympy import *
import csv
import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import minimize
from scipy.interpolate import interp1d
from Qscheme.simulate import array_eins_from_params_product

class Fit_Base:
    def __init__(self,params_marked):
        self.params_marked = params_marked
        self.mut_params_x0_map = [i for i,x_marked in enumerate(params_marked) if x_marked[1] == 1]
        self.mut_params_x0 = [params_marked[i][0] for i in self.mut_params_x0_map]
        self.fixed_params_map = [i for i,x_marked in enumerate(params_marked) if x_marked[1] == 0]
        self.fixed_params = [params_marked[i][0] for i in self.fixed_params_map]

        self.data_dir = "C:\\Users\\user\\Documents\\GitHub\\Qdev\\qubit_simulations\\Qutip&SymPy simulations of 3 junction Cshunted flux qbit based on 3\\Design v.1.0\\"
        self.data_files_names = ["E1E0WPD.csv","E2E0halfWPD.csv","E2E1WPD.csv"]
        self.data_names = [x.split('.')[0] for x in self.data_files_names]
        self.eigEngs_current_step_dict = {}
        self.eigEngs_data_dict = {}
        self.eigEngs_data_splineFits_dict = {}
        self.eigEngs_data_inter_dict = {}
    
        self.visualize = True
        
    def import_data(self):
        for data_name, data_file_name in zip(self.data_names,self.data_files_names):
            with open(self.data_dir + data_file_name, "r" ) as file:
                rows = list(csv.reader(file))
                self.eigEngs_data_dict[data_name] = np.array( rows, dtype=np.float64 ).T
                x,y = self.eigEngs_data_dict[data_name]
                self.eigEngs_data_splineFits_dict[data_name] = interp1d(x,y,kind="cubic")
        
        lb_list = [min(self.eigEngs_data_dict[data_name][0]) for data_name in self.data_names]
        rb_list = [max(self.eigEngs_data_dict[data_name][0]) for data_name in self.data_names]
        
        self.spline_data_lb = max(lb_list)
        self.spline_data_rb = min(rb_list)
                
            
    def run_fit(self,pts_N):
        fs = np.linspace(self.spline_data_lb,self.spline_data_rb, pts_N)
        self._calculate_splines_for_new_fs(fs)
        
        if( self.visualize == True ):
            plt.ion()
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111)
            self.data_lines = []
            self.sim_lines = []
            for i,data_name in enumerate(self.data_names):
                Ei = self.eigEngs_data_inter_dict[data_name]
                self.data_lines.append(self.ax.plot(self.fs,Ei,label=data_name)[0])
                self.sim_lines.append(self.ax.plot([],[])[0])
            self.ax.legend()
            self.ax.grid()
            self.func_to_minimize_wrap(self.mut_params_x0)
            self.iter_callback(self.mut_params_x0)
            
        
        tol = 1e-6
        if( self.visualize == True ):
            callback_method = self.iter_callback
        else:
            callback_method = None
        
        self.fit_result = minimize(self.func_to_minimize_wrap, self.mut_params_x0,
                        method='Powell', tol=None, callback=callback_method, 
                        options={'disp': False, 'return_all': False, 'maxiter': None,
                                 'direc': None, 'maxfev': None,
                                 'xtol': tol, 'ftol': tol}
                       )
        self.iter_callback(self.fit_result.x)
        return self.fit_result
    
    def _calculate_splines_for_new_fs(self,fs):
        self.fs = fs
        for data_name in self.data_names:
            self.eigEngs_data_inter_dict[data_name] = np.array(
                    self.eigEngs_data_splineFits_dict[data_name](self.fs)
                    )
    
    def func_to_minimize_wrap(self,mut_params):
        params = [None for i in self.params_marked]
        for i,i_map in enumerate(self.mut_params_x0_map):
            params[i_map] = mut_params[i]
        for i,i_map in enumerate(self.fixed_params_map):
            params[i_map] = self.fixed_params[i]
                    
        return self.func_to_minimize(params)
        
    def func_to_minimize(self,params):
        raise NotImplementedError
        
    
    def iter_callback(self,xk):
        print("callback ", xk)
        
        for i,data_name in enumerate(self.data_names):
            self.sim_lines[i].set_xdata(self.fs)
            self.sim_lines[i].set_ydata(self.eigEngs_current_step_dict[data_name])
        self.ax.relim()
        self.ax.autoscale_view()
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        
        

class Fit_Flux_Csh_3JJ(Fit_Base):
    def func_to_minimize(self, params):
        x_alpha = params[0]
        x_beta = params[1]
        E_C = params[2]
        E_J = params[3]
        Csh = params[4]
        alpha = params[5]
        cooperN = params[6]

        Evs, Engs = array_eins_from_params_product(x_alpha*(self.fs - x_beta),E_C,E_J,Csh,alpha,cooperN,4)
        
        self.eigEngs_current_step_dict[self.data_names[0]] = Engs[:,1] - Engs[:,0]
        self.eigEngs_current_step_dict[self.data_names[1]] = (Engs[:,2] - Engs[:,0])/2
        self.eigEngs_current_step_dict[self.data_names[2]]= Engs[:,2] - Engs[:,1]
        
        
        d = 0
        for data_name in self.data_names:
            E_data = self.eigEngs_data_inter_dict[data_name]
            E_sim = self.eigEngs_current_step_dict[data_name]
            d += np.sum( (E_data - E_sim)**2 )
        print(1e5*d)
        return 1e5*d
    

if( __name__ == "__main__" ):
    alpha = 0.43 # ratio of the smaller and larger josephson junctions areas
    E_C = 8.56735249e+00  # GHz, roughly corresponding to C = 4.5 fF
    E_J = 1.22570055e+02 # GHz
    Csh = 10.0 # in units of C (C equals to larger josephson junctions capacitance)
    f = 0.5 # flux
    cooper_N = 15
    x_alpha =  4.84244299e-02
    x_beta = -1.63358728e+01
    pts_N = 25
    
    params_marked = [[x_alpha,1], [x_beta,1], [E_C,1], [E_J,1],
              [Csh,0], [alpha,0], [cooper_N,0]]
    
    fit = Fit_Flux_Csh_3JJ(params_marked)
    fit.import_data()
    print(fit.run_fit(pts_N))
    
    
