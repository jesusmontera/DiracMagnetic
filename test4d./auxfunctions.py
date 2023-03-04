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

    
def make2DGaussian(N,L,k,pos,sigma=0.06):
    # p is momentum, k is the wave number that is the# spatial frequency
    # with respect to spatial extent of the simulation
    p = k * 2.0 * np.pi / L        
    X = L*np.linspace(-0.5, 0.5 - 1.0/N, N)
    X1, Y1 = np.meshgrid(X, X)
    pf = np.exp(1j*p[0]*X1 + 1j*p[1]*Y1)
    gaussian = pf * np.exp(-((X1/L+pos[0]/2.)/sigma)**2/2.0                                    
                  - ((Y1/L-pos[1]/2.)/sigma)**2/2.0 )
    #return gaussian 
    return gaussian/np.sqrt(np.sum(gaussian*np.conj(gaussian)))
    

def getEnergyEigenSpinors(N,L,k,m=1.):
    
    C = 137.036 # Speed of light
    #C=1.0
    mc = m*C
    
    px, py, pz = 2.0*k[0]*np.pi/L, 2.0*k[1]*np.pi/L, 2.0*k[2]*np.pi/L
    
    
    p2 = px**2 + py**2 + pz**2
    p = np.sqrt(p2)
    omega = np.sqrt(mc*mc + p2)

    den1 = p*np.sqrt((mc - omega)**2 + p2)
    den2 = p*np.sqrt((mc + omega)**2 + p2) # from 2d
        
    #zeros = np.zeros(N_DIM*[N], dtype=np.complex128)
    omega = np.sqrt(mc*mc + p2) # Corresponds to E/c

    # Temporary variable used for
    # some denominators in the negative energy solutions
    den1 = p*np.sqrt((mc - omega)**2 + p2)
    # Used for denominators in the positive energy solutions
    den2 = p*np.sqrt((mc + omega)**2 + p2)
    # energy eigenstates in bra form
    # Negative energy solutions
    neg_eig1 = [pz*(mc - omega)/den1, 
                (mc*px - 1.0j*mc*py - (px - 1.0j*py)*omega)/den1,
                p2/den1, 0.]
    neg_eig2 = [(mc*px + 1.0j*mc*py + (-px - 1.0j*py)*omega)/den1,
                -pz*(mc - omega)/den1, 0., p2/den1]

    # Positive energy solutions
    
    pos_eig1 = [pz*(mc + omega)/den2, 
                (mc*px - 1.0j*mc*py + (px - 1.0j*py)*omega)/den2,
                p2/den2, 0.]
    pos_eig2 = [(mc*px + 1.0j*mc*py + (px + 1.0j*py)*omega)/den2,
                -pz*(mc + omega)/den2, 0., p2/den2]
        
    return pos_eig1,pos_eig2, neg_eig1,neg_eig2
def initspinors3D(pos_eig1,pos_eig2,neg_eig1,neg_eig2, spin=None, bPositive=True):    

    if spin is None:
        # if no initial spin is requiered  then spin will point to momentum direction
        spinorsketpos = np.array([pos_eig1[1], 0.0, 0.0, pos_eig1[2]])
        sinfo =" (spin coincides with Kinetic)"
    else:        
        # if initial spin is requiered  make spinors for that specific spin
        # convert spinors to ket form to merge spin to a momentum
        sinfo =" from specific spin" + str(spin) + "\n\t"
        if bPositive:
                            
            sinfo =" positive eigen energies "
            c1 = pos_eig1 @ np.array([spin[0], spin[1], 0.0, 0.0]) # inner product
            c2 = pos_eig2 @ np.array([spin[0], spin[1], 0.0, 0.0])
            spinorsketpos = c1 * np.conj(pos_eig1) + c2 * np.conj(pos_eig2)
        else:
            sinfo =" negative eigen energies "
            c1 = neg_eig1 @ np.array([ 0.0, 0.0,spin[0], spin[1]]) # inner product
            c2 = neg_eig2 @ np.array([ 0.0, 0.0,spin[0], spin[1]])
            spinorsketpos = c1 * np.conj(neg_eig1) + c2 * np.conj(neg_eig2)

    
    spinorsketpos = spinorsketpos / np.linalg.norm(spinorsketpos)        
    print(sinfo)            
    init_spinor = [spinorsketpos[0],spinorsketpos[1],spinorsketpos[2],spinorsketpos[3] ]
    return init_spinor
def makeGaussian2P2D(N,L,pos1,pos2,k1,k2):
    X = L*np.linspace(-0.5, 0.5 - 1.0/N, N)
    Y1, Y2, X1, X2 = np.meshgrid(X, X, X, X)
    k_to_p = (2.0 * np.pi / L)
    px1 = k1[0] * k_to_p
    px2 = k2[0] * k_to_p
    py1 = k1[1] * k_to_p
    py2 = k2[1] * k_to_p
    pf = np.exp(1j*px1*X1 + 1j*py1*Y1 + 1j*px2*X2 + 1j*py2*Y2)
    SIGMA = np.complex128(0.056568)
    xa0=-pos1[0]/2.
    xb0=pos2[0]/2.
    print("xa0",xa0)
    print("xb0",xb0)
    wavefunc = pf * np.exp(-((X1/L+xa0)/SIGMA)**2/2.0
                  - ((X2/L-xb0)/SIGMA)**2/2.0
                  - ((Y1/L)/SIGMA)**2/2.0
                  - ((Y2/L)/SIGMA)**2/2.0
                  )
    wavefunc = wavefunc/np.sqrt(np.sum(wavefunc*np.conj(wavefunc)))
    extent=[Y1[0, 0, 0, 0], Y1[0, -1, 0, 0],
            X1[0, 0, 0, 0], X1[0, 0, -1, 0]]
    print(extent)
    return (wavefunc - np.transpose(wavefunc, (1, 0, 3, 2)))/np.sqrt(2.0), extent
    
def makeGaussian2P2Dold(N,L,pos1,pos2,k1,k2):    

    k_to_p = -(2.0 * np.pi / L)
    px1 = k1[0] * k_to_p
    px2 = k2[0] * k_to_p
    py1 = k1[1] * k_to_p
    py2 = k2[1] * k_to_p
    σ = np.complex128(0.056568)
    
    
    x1, x2,y1,y2 = np.mgrid[ -L/2: L/2:N*1j, -L/2: L/2:N*1j, -L/2: L/2:N*1j, -L/2: L/2:N*1j]

    pf = np.exp(1j*px1*x1 + 1j*py1*y1 + 1j*py2*x2 + 1j*px2*y2)
    
    wf= pf * np.exp( -1/(4* σ**2) * ((x1-pos1[0])**2+(y1-pos1[1])**2+
                                     (x2-pos2[0])**2 +(y2-pos2[1])**2)) / np.sqrt(2*np.pi* σ**2)
    extent=[-1,1,-1,1]

    
    return  wf,extent

    
    
    
