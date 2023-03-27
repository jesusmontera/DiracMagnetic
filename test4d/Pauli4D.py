import numpy as np

from splitstep import SplitStepMethod 
import scipy.constants as const
from auxfunctions import pauli4x4Matrixs
from potentials import applyrotmatrix4D #, makebarrier4d,makeBpotential4D
from numba import jit


class Pauli4D():
    
    def __init__(self):                
        self.M = 1.0 # Mass of the particle
        self.C =  137.036 # Speed of light
        self.HBAR = 1.0
        self.N=32
        self.L=2.
        self.dt=0.001
        self.U = None                
        self.psi=None
        self.V=None 
        self.numframes=0
        self.psi2d=None      
        self.rotmatrix=None
        self.exp_magnetic=None # B magnetic potential
        self.pauli4x, self.pauli4y, self.pauli4z = pauli4x4Matrixs()
    def clear(self):
        
        self.psi=None
        self.U=None
        self.exp_magnetic=None
        self.numframes=0
    def makerotMatrix(self):
        # this 4d rotation spinor matrix simulates a magnetic field
        N= self.N
        rotmatrix=np.zeros([N,N,N,N,4,4],dtype=np.complex128)
    
        for x in range(N):            
            for y in range(N):                
##                if x<N//2:
##                    σ = 0.0
##                else:
##                    σ = np.pi*y*0.7 /N
                σ = np.pi*y*0.7 /N  # coment this line and uncoment the upper to apply
                                    # the magnetic field(rotations) to only one particle
                    
                c=np.cos(σ/2)
                s= np.sin(σ/2)
                r2= np.array([[c,s],[-s,c]])
                r4 = np.kron(r2,r2) # this convert spinors rotations from 2d to 4d
                # the next 2 loops loops is to convert the rotations from 2d  to 4d
                for a in range(N):
                    for b in range(N):                
                        rotmatrix[x][y][a][b]=r4
                
        return rotmatrix
    
    
       
    def makewavefuncion(self,N,L,pos,k):
        X = L*np.linspace(-0.5, 0.5 - 1.0/N, N)
        X1, X2, Y1, Y2 = np.meshgrid(X, X, X, X)

        SIGMA = np.complex128(0.056568)
                
        p = k * ( 2.0 * np.pi / L )
        pf = np.exp(1j*p[0]*X1 + 1j*p[1]*Y1 + 1j*p[2]*X2 + 1j*p[3]*Y2)

        pos /=2.
        x1=-pos[0]
        x2=pos[2]
        y1=-pos[1]
        y2=pos[3]
        wf = pf*np.exp(-((X1/L+x1)/SIGMA)**2/2.0
                          - ((X2/L-x2)/SIGMA)**2/2.0
                          - ((Y1/L+y1)/SIGMA)**2/2.0
                          - ((Y2/L-y2)/SIGMA)**2/2.0
                          )
        wf = wf/np.sqrt(np.sum(wf*np.conj(wf)))

        wf = (wf - np.transpose(wf, (1, 0, 3, 2)))/np.sqrt(2.0)
        return wf
                
        
    def initialize(self,N,L,DT,pos,k,initial_spin):

        
        self.N  = N
        self.L=L
        self.dt=DT        

        np.seterr(under="ignore")                
        self.numframes=0
        self.psi2d= np.zeros([2,4,N,N], dtype=np.complex128)
        
                    
        wavefunc = self.makewavefuncion(N,L,pos,k)
        
                       
        ones = np.ones([N,N,N,N], dtype=np.complex128)            
        self.psi=wavefunc * np.multiply.outer(initial_spin, ones)
        
##        n = np.linalg.norm(self.psi)
##        if n != 0:
##            self.psi /= n
            
        
        Lm=L*5.29177210903E-11
        DTm=DT*2.4188843265857E-17
        if self.V is None:
            #self.V= makebarrier4d(N)            
            self.V=np.zeros([N,N,N,N])
        self.rotmatrix = self.makerotMatrix()                
        
        self.U = SplitStepMethod(self.V, (Lm, Lm,Lm, Lm), DTm)

        
    def getProb(self):
        return sum([ np.einsum('ijkl->jk', np.abs(self.psi[i]))for i in range(4)])    

    
    def dostep(self):
        
        for i in range(4):
            self.psi[i] = self.U(self.psi[i])
        if self.rotmatrix is not None:            
            applyrotmatrix4D(self.N, self.psi, self.rotmatrix)
            
        if self.exp_magnetic is not None:            
            self.psi = np.einsum('ij...,j...->i...', self.exp_magnetic, self.psi)
        self.numframes+=1
    def psi4dto2d(self):
        for i in range(4):
            self.psi2d[0][i] = np.einsum('jqik->ij', self.psi[i])# /self.sumpsi2d[1][i]         
            self.psi2d[1][i]= np.einsum('kiqj->ij', self.psi[i]) #/ self.sumpsi2d[0][i]
    
    def getSpinExpecValue(self):
        self.psi4dto2d()    
        sex=np.zeros([2,3],dtype=np.complex128)
        
        for i in range(2):
            sa = np.sum(np.abs(self.psi2d[i])**2)
            cp=np.conj(self.psi2d[i])
            sex[i][0] = np.sum(cp*np.einsum('ij,j...->i...', self.pauli4x, self.psi2d[i]))/sa
            sex[i][1] = np.sum(cp*np.einsum('ij,j...->i...', self.pauli4y, self.psi2d[i]))/sa
            sex[i][2] = np.sum(cp*np.einsum('ij,j...->i...', self.pauli4z, self.psi2d[i]))/sa
            #sex[i] /= np.linalg.norm(sex[i])
    
        return np.real(sex)
    
        

        
            
        
    
        
