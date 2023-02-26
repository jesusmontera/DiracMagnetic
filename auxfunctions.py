import numpy as np
import math
from numba import jit
##########################################################################
#    auxiliary function input qubit spin output it's bloch vector
##########################################################################
def getBlochVector(spin):
    #compute density matrix
    rho=np.zeros([2,2],np.complex128)
    for i in range(2):
        for k in range(2):
            if i == k:
                rho[i][k]=spin[i]**2
            else:
                rho[i][k]=spin[i] * np.conj(spin[k])
            
    a = rho[0, 0]
    b = rho[1, 0]
    x = 2.0 * b.real
    y = 2.0 * b.imag
    z = 2.0 * abs(a) - 1.0
    bloch=np.array([x,y,z])        
    return bloch
def pauli4x4Matrixs():    
    p2x = np.array([[0.0, 1.0], [1.0, 0.0]],np.complex128)
    p2y = np.array([[0,-1j],[1j,0]],np.complex128)
    p2z = np.array([[1,0],[0,-1]],np.complex128)

    p4x = np.block([[p2x, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2x]])
    p4y = np.block([[p2y, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2y]])
    p4z = np.block([[p2z, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2z]])
    
    return p4x,p4y,p4z
        

@jit(nopython=True)
def spinDotBdirac(N: int, DT: float, wf: np.ndarray,B: np.ndarray, spinbloch: np.ndarray):    

    for x in range(N):
        for y in range(N):
            for z in range(N):
                for i in range(4):
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
    #p = k * 2.0 * np.pi / L
    p = k * np.pi / L
    print("gaussian momentum",p)
    X, Y , Z  = np.mgrid[ -L/2: L/2:N*1j, -L/2: L/2:N*1j, -L/2: L/2:N*1j]
    
    pf= np.exp(1j*p[0]*X + 1j*p[1]*Y +1j*p[2]*Z )
    
    gaussian = pf * np.exp( -1./(4* sigma**2) * ( (X-pos[0])**2 + (Y-pos[1])**2  + (Z-pos[2])**2 ))  * np.sqrt(2*np.pi* sigma**2) 
    return gaussian #gaussian/np.sqrt(np.sum(gaussian*np.conj(gaussian)))

    
    
    
