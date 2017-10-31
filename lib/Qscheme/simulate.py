# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 20:17:14 2017

@author: user
"""

import numpy as np
import qutip as qp
import itertools

e = 1.6e-19 # electron charge SI
hbar = 1.055e-34 # plank constant SI (this is \hbar)
h = hbar*(2*np.pi)

alpha = 0.43 # ratio of the smaller and larger josephson junctions areas
E_C = 17.6 # GHz, roughly corresponding to C = 4.5 fF
C = 4*e**2/2/(2*np.pi*E_C*hbar*10**9)
E_J = 36.2/alpha # GHz
Csh = 51e-15/C # in units of C (C equals to larger josephson junctions capacitance)

Z_r = 50 # Ohm - resonator line impendance
freq_r = 6.00 # GHz resonator frequency
omerga_r = freq_r*(2*np.pi)
beta = 0.043 # coupling parameter
#rms single-photon voltage in resonator
V_r0 = freq_r*np.sqrt(2*h*Z_r)

def get_inv_dl_caps_matrix( Csh, alpha ):
    return 1/(2*(alpha + Csh) + 1)*np.array( [[1+alpha+Csh,alpha+Csh],[alpha + Csh,1 + alpha + Csh]] )

def H_c( Ec, Csh, alpha, N ):
    islands_N = 2
    charge_ops = qp.tensor(qp.charge(N),
                           qp.identity(2*N+1)), \
                 qp.tensor(qp.identity(2*N+1),
                           qp.charge(N))
    
    H_c = qp.tensor(qp.qzero(2*N+1), qp.qzero(2*N+1))
    inv_dl_caps_matrix = get_inv_dl_caps_matrix( Csh, alpha )
    for i in range(0,islands_N):
        for j in range(0,islands_N):
            H_c += Ec*charge_ops[i]*inv_dl_caps_matrix[i,j]*charge_ops[j]
    return H_c

def raising(N):
    return qp.Qobj(np.diag(np.ones(2*N),-1))
def lowering(N):
    return qp.Qobj(np.diag(np.ones(2*N),1))
# tunneling operator is available in QuTip and his result equals to "rising(N) + lowering(N)"

def H_J( E_J, alpha, phase_ext, N ):
    H1 = 0.5*( np.exp(1j*phase_ext)*qp.tensor(lowering(N),raising(N)) \
                         + np.exp(-1j*phase_ext)*qp.tensor(raising(N),lowering(N)) )
    H2 = qp.tensor(qp.identity(2*N+1), 0.5*qp.tunneling(2*N+1,1))
    H3 = qp.tensor(0.5*qp.tunneling(2*N+1,1), qp.identity(2*N+1))
        
    tunneling_ops = (H1,H2,H3)                  
    H_J = qp.tensor(qp.qzero(2*N+1), qp.qzero(2*N+1))
    E_Js = (-1)*np.array([alpha*E_J,E_J,E_J])
    for idx, tunneling_op in enumerate(tunneling_ops):
        H_J += E_Js[idx]*tunneling_op
    return H_J

def array_eins(f_0,f_1,Nf,E_C,E_J,Csh,alpha,cooper_N,eigvals_N):
    arr_evals = []
    arr_einvects = []
    fs = np.linspace(f_0,f_1,Nf)
    for f in fs: 
        print('\r{:g} '.format(f), end='')
        H = H_J(E_J,alpha,f,cooper_N) + H_c(E_C,Csh,alpha,cooper_N)
        evals, einvects = H.eigenstates(sparse=True,eigvals=eigvals_N)
        arr_evals.append(evals)
        arr_einvects.append(einvects)
    return np.array(arr_einvects), np.array(arr_evals), np.array(fs)

def plot_eigenergies(engs,fs):
    for idx in range(1,engs.shape[1]):
        plt.plot(fs,engs[:,idx]-engs[:,0], label=r'$E_{%s}-E_{%s}$'%(idx,0))
    plt.plot( fs,engs[:,2]-engs[:,1], label=r"$E_2 - E_1$" )
    plt.plot( fs,(engs[:,3]-engs[:,0])/2, label=r"$(E_3 - E_0)/2$" )
    plt.plot( fs,(engs[:,2]-engs[:,0])/2, label=r"$(E_2 - E_0)/2$" )
    plt.legend()
    plt.grid()
    
def plot_eigenergies_custom(engs,fs):
    for idx in range(1,3):
        plt.plot(fs,engs[:,idx]-engs[:,0], label=r'$E_{%s}-E_{%s}$'%(idx,0))
    plt.plot( fs,engs[:,2]-engs[:,1], label=r"$E_2 - E_1$" )
    plt.plot( fs,(engs[:,3]-engs[:,0])/2, label=r"$(E_3 - E_0)/2$" )
    plt.plot( fs,(engs[:,2]-engs[:,0])/2, label=r"$(E_2 - E_0)/2$" )
    plt.legend()
    plt.grid()


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
        
    return np.array(arr_einvects), np.array(arr_evals)

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
        
    return np.array(arr_einvects), np.array(arr_evals)

def array_coupling_from_parameters_product(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts):
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

def array_coupling_from_parameters_lists(fl_pts,E_C_pts,E_J_pts,Csh_pts,alpha_pts,cooper_N_pts):
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