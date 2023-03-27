import numpy as np
import math
from numba import jit

@jit(nopython=True)    
def applyDotMatrix3D( N: int , psi: np.ndarray, Bpot : np.ndarray,bv : np.ndarray):    
    sv= np.zeros(2,dtype=np.complex128)
    for x in range (N):
        for y in range (N):
            for z in range(N):                
                sv[0]= psi[0][x][y][z]
                sv[1]= psi[1][x][y][z]
                
##                bv[0][0]=1+1j
##                bv[0][1]=1+2j
##                bv[1][0]=1+3j
##                bv[1][1]=1+5j

                bv[0][0]=Bpot[0][0][x][y][z]
                bv[0][1]=Bpot[0][1][x][y][z]
                bv[1][0]=Bpot[1][0][x][y][z]
                bv[1][1]=Bpot[1][1][x][y][z]

                
                sv = np.dot(sv,bv)
                psi[0][x][y][z]=sv[0]
                psi[1][x][y][z]=sv[1]
                
                


##########################################################################
#   construct potential from B magnetic vector's matrix 
##########################################################################

def makeBpotential(dt, B):
    # aplyed with
    #wave_func = np.einsum('ij...,j...->i...', exp_magnetic, wave_func).
    if B is None:
        return None    
    q = 1. #1.602176634e-19 #Coulombs
    a = q*np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)*dt
    d=np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    d[d == 0] = 1e-40    
    
    nx, ny, nz = B/d
    exp_magnetic = np.array([[np.cos(a) + 1.0j*nz*np.sin(a),
                              1.0j*(nx - 1.0j*ny)*np.sin(a)],
                             [1.0j*(nx + 1.0j*ny)*np.sin(a),
                              np.cos(a) - 1.0j*nz*np.sin(a)]])
    print("exp_magnetic",np.shape(exp_magnetic))
    return exp_magnetic

##########################################################################
#    input qubit spin output it's bloch vector
##########################################################################
                
def getBlochVector(spin):
    
    rho = np.outer(spin, np.conj(spin)) #compute density matrix
            
    a = rho[0, 0]
    b = rho[1, 0]
    x = 2.0 * b.real
    y = 2.0 * b.imag
    z = 2.0 * abs(a) - 1.0
    bloch=np.array([x,y,z])        
    return bloch
def pauli2x2Matrixs():    
    p2x = np.array([[0.0, 1.0], [1.0, 0.0]],np.complex128)
    p2y = np.array([[0,-1j],[1j,0]],np.complex128)
    p2z = np.array([[1,0],[0,-1]],np.complex128)
    return p2x,p2y,p2z
def pauli4x4Matrixs():    
    p2x,p2y,p2z =  pauli2x2Matrixs()
    
    p4x = np.block([[p2x, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2x]])
    p4y = np.block([[p2y, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2y]])
    p4z = np.block([[p2z, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2z]])
    
    return p4x,p4y,p4z
        

##########################################################################
#   apply B vector's matrix to each of the four Dirac psi's  or Paulis two psi's
##########################################################################

@jit(nopython=True)
def spinDotB(N: int, DT: float, wf: np.ndarray,B: np.ndarray, spinbloch: np.ndarray):    

    ni=len(wf) # dim's psi (2 pauli) ( 4 dirac)
    
    for x in range(N):
        for y in range(N):
            for z in range(N):
                for i in range(ni):
                    R = wf[i][x][y][z].real
                    I = wf[i][x][y][z].imag                                            
                    # imag                    
                    dd = spinbloch.dot(B[x][y][z]) #* prob #(/maxprob)
                    Iinc =  dd * R * DT 
                    #real
                    Rdec = dd * I * DT
                    inc= -Rdec + Iinc *1j
                    #apply B contribution
                    wf[i][x][y][z] += inc    
    
def make3DGaussian(N,L,k,pos,sigma=0.06):
    # p is momentum, k is the wave number that is the# spatial frequency
    # with respect to spatial extent of the simulation
    p = k * 2.0 * np.pi / L
    #p = k * np.pi / L
    print("gaussian momentum",p)
    X, Y , Z  = np.mgrid[ -L/2: L/2:N*1j, -L/2: L/2:N*1j, -L/2: L/2:N*1j]
    
    pf= np.exp(1j*p[0]*X + 1j*p[1]*Y +1j*p[2]*Z )
    
    gaussian = pf * np.exp( -1./(4* sigma**2) * ( (X-pos[0])**2 + (Y-pos[1])**2  + (Z-pos[2])**2 ))  * np.sqrt(2*np.pi* sigma**2) 
    return gaussian #gaussian/np.sqrt(np.sum(gaussian*np.conj(gaussian)))

    
    
    
