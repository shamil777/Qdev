import qutip as qp
import numpy as np

def charge_ops( nodes_N,cooper_N ):
    charge_ops = []
    for i in range(nodes_N):
        to_tensor_list = []
        for j in range(nodes_N):
            if( j != i ):
                to_tensor_list.append( qp.identity(2*cooper_N+1) )
            else:
                to_tensor_list.append( qp.charge(cooper_N) )
        charge_ops.append( qp.tensor(*to_tensor_list) )
    return charge_ops

# tunneling operator is available in QuTip and his result equals to "rising(N) + lowering(N)"
def raising(cooper_N):
    return qp.Qobj(np.diag(np.ones(2*cooper_N),-1))
    
def lowering(N):
    return qp.Qobj(np.diag(np.ones(2*N),1))

def raising_lowering_ops(nodes_N,cooper_N):
    raising_ops = []
    lowering_ops = []
    for i in range(nodes_N):
        to_tensor_list_r = []
        to_tensor_list_l = []
        for j in range(nodes_N):
            if( j != i ):
                to_tensor_list_r.append( qp.identity(2*cooper_N+1) )
                to_tensor_list_l.append( qp.identity(2*cooper_N+1) )
            else:
                to_tensor_list_r.append( raising(cooper_N) )
                to_tensor_list_l.append( lowering(cooper_N) )
        raising_ops.append( qp.tensor(*to_tensor_list_r) )
        lowering_ops.append( qp.tensor(*to_tensor_list_l) )
        
    return raising_ops,lowering_ops