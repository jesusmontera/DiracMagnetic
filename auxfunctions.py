import numpy as np
import math
from numba import jit
##########################################################################
#    auxiliary function input qubit spin output it's bloch vector
##########################################################################
def getBlochVector(statevector):
    alpha_real = statevector[0].real
    alpha_imag = statevector[0].imag

    if alpha_real == 0: alpha_real = 1e-20
    
    alpha_theta = math.atan(alpha_imag/alpha_real) 

    if alpha_real < 0 and alpha_imag > 0 :
        alpha_theta = math.atan(alpha_imag/alpha_real) + math.pi

    if alpha_real < 0 and alpha_imag < 0 :
        alpha_theta = math.atan(alpha_imag/alpha_real) + math.pi
    r_alpha = math.sqrt((alpha_real**2) + (alpha_imag**2))

    beta_real = statevector[1].real
    beta_imag = statevector[1].imag
    if beta_real == 0: beta_real=1e-20
    beta_theta = math.atan(beta_imag/beta_real) 

    if beta_real < 0 and beta_imag > 0 :
        beta_theta = math.atan(beta_imag/beta_real) + math.pi

    if beta_real < 0 and beta_imag < 0 :
        beta_theta = math.atan(beta_imag/beta_real) + math.pi
    r_beta = math.sqrt((beta_real**2) + (beta_imag**2))        
    phi = beta_theta - alpha_theta
    theta = 2 * np.arccos(r_alpha)
    blochvector= np.zeros(3)
    blochvector[0] = math.sin(theta)*math.cos(phi)
    blochvector[1] = math.sin(theta)*math.sin(phi)
    blochvector[2] = math.cos(theta)
    return blochvector
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

    
    
    
