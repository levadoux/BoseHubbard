import math
import numpy as np
from scipy.sparse.linalg import eigsh

# LABELING FUNCTION

def size_Fock_space(N,M):
    return math.comb(N+M-1,N)

def state_to_positions(N,M,state):
    positions = [0]*N
    start = 0
    for j in range(M):
        to_placed = state[M-1-j]
        positions[start:start+to_placed]=[M-j]*to_placed
        start += to_placed
    return positions

def positions_to_label(N,M,positions):
    label = 1
    for j in range(N):
        label += size_Fock_space(j+1,M-positions[j])
    return label

def labeling(N,M,state):
    label = positions_to_label(N,M,state_to_positions(N,M,state))
    return label

# UNLABELING FUNCTION

def label_to_positions(N,M,label):
    positions = [0]*N
    rest = label
    for j in range(N):
        exp = 0
        while size_Fock_space(N-j,exp) < rest:
            exp += 1
        m = M+1-exp # exp = M-m+1 at the end of the while loop
        positions[N-1-j] = m
        rest -= size_Fock_space(N-j,M-m)
    return positions

def positions_to_state(N,M,positions):
    state = [0]*M
    for site in range(M):
        state[site] = positions.count(site+1) # sites are indexed from 1 to M
    return state

def unlabeling(N,M,label):
    state = positions_to_state(N,M,label_to_positions(N,M,label))
    return state

# CREATION_ANNIHILATION FUNCTION (defined on Fock states)

def creation_annihilation(i,j,state):
    new_state = state.copy()
    n_i = new_state[i-1]
    n_j = new_state[j-1]
    if i==j:
        coeff = state[i-1]
    elif n_j>0:
        new_state[i-1] = n_i+1
        new_state[j-1] = n_j-1
        coeff = math.sqrt(n_j)*math.sqrt(n_i+1)
    else:
        coeff = 0
    return new_state, coeff

# BOSE-HUBBARD HAMILTONIAN

def BoseHubbard_Hj(N,M):
    dim = size_Fock_space(N,M)
    Hj = np.zeros((dim,dim))
    for label in range(1,dim+1):
        state = unlabeling(N,M,label)
        for j in range(1,M):
            new_state1, coeff1 = creation_annihilation(j,j+1,state)
            new_state2, coeff2 = creation_annihilation(j+1,j,state)
            if coeff1 != 0:
                new_label1 = labeling(N,M,new_state1)
                Hj[new_label1-1,label-1] = coeff1
            if coeff2 != 0:
                new_label2 = labeling(N,M,new_state2)
                Hj[new_label2-1,label-1] = coeff2
        if M>2:
            new_state1, coeff1 = creation_annihilation(M,1,state)
            new_state2, coeff2 = creation_annihilation(1,M,state)
            if coeff1 != 0:
                new_label1 = labeling(N,M,new_state1)
                Hj[new_label1-1,label-1] = coeff1
            if coeff2 != 0:
                new_label2 = labeling(N,M,new_state2)
                Hj[new_label2-1,label-1] = coeff2
    return Hj

def BoseHubbard_Hu(N,M):
    dim = size_Fock_space(N,M)
    Hu = np.zeros((dim,dim))
    for label in range(1,dim+1):
        state = unlabeling(N,M,label)
        for n_s in state:
            Hu[label-1,label-1] += 0.5*n_s*(n_s-1)
    return Hu

# CONSTRUCTION OF THE ONE-BODY REDUCED DENSITY MATRIX

def Single_Hopping(N,M,i,j):
    dim = size_Fock_space(N,M)
    H_ij = np.zeros((dim,dim))
    for label in range(1,dim+1):
        state = unlabeling(N,M,label)
        new_state, coeff = creation_annihilation(i,j,state)
        new_label = labeling(N,M,new_state)
        H_ij[new_label-1,label-1] = coeff
    return H_ij


def OBDM_ground(N,M,U,J):
    rho1 = np.zeros((M,M))
    H = -J*BoseHubbard_Hj(N,M)+U*BoseHubbard_Hu(N,M)
    ground = eigsh(H,1,which='SA')[1]
    for i in range(1,M+1):
        for j in range(1,M+1):
            H_ij = Single_Hopping(N,M,i,j)
            rho1[i-1,j-1] = (ground.T @ H_ij @ ground)[0,0]
    return rho1

# CONSTRUCTION OF THE TWO-BODY REDUCED DENSITY MATRIX

def crea_crea_ann_ann(k,l,j,i,state):
    new_state = state.copy()
    n_i = new_state[i-1]
    n_j = new_state[j-1]
    if i==j and n_i>1:
        new_state[i-1] = n_i-2
        n_l = new_state[l-1]
        new_state[l-1] +=1
        n_k = new_state[k-1]
        new_state[k-1] +=1
        coeff = math.sqrt(n_i)*math.sqrt(n_i-1)*math.sqrt(n_l+1)*math.sqrt(n_k+1)
    elif i!=j and n_i>0 and n_j>0:
        new_state[i-1] = n_i-1
        new_state[j-1] = n_j-1
        n_l = new_state[l-1]
        new_state[l-1] +=1
        n_k = new_state[k-1]
        new_state[k-1] +=1
        coeff = math.sqrt(n_i)*math.sqrt(n_j)*math.sqrt(n_l+1)*math.sqrt(n_k+1)
    else:
        coeff = 0
    return new_state, coeff

def Pair_Hopping(k,l,j,i,N,M):
    dim = size_Fock_space(N,M)
    Hijkl = np.zeros((dim,dim))
    for label in range(1,dim+1):
        state = unlabeling(N,M,label)
        new_state, coeff = crea_crea_ann_ann(k,l,j,i,state)
        new_label = labeling(N,M,new_state)
        Hijkl[new_label-1,label-1] = coeff
    return Hijkl

def TBDM_ground(N,M,U,J):
    rho2 = np.zeros((M*M,M*M))
    H = -J*BoseHubbard_Hj(N,M)+U*BoseHubbard_Hu(N,M)
    ground = eigsh(H,1,which='SA')[1]
    for i in range(1,M+1):
        for j in range(1,M+1):
            for k in range(1,M+1):
                for l in range(1,M+1):
                    Hijkl = Pair_Hopping(k,l,j,i,N,M)
                    rho2[(i-1)*M+(j-1),(k-1)*M+(l-1)]=(ground.T @ Hijkl @ ground)[0,0]
    return rho2







