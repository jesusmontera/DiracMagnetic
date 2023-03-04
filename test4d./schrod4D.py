import numpy as np
from splitstep import SplitStepMethod
from auxfunctions import make2DGaussian


class schrod4D():
    
    def __init__(self,N,V):                
        self.M = 1.0 # Mass of the particle
        self.C =  137.036 # Speed of light
        self.HBAR = 1.0
        self.N=N
        self.dt=0.001
        self.U = None        
        self.psi4d=None
        self.psi2dA=None
        self.psi2dB=None
        self.V = V
        self.sumA=0
        self.sumB=0
    def initialize(self,L,DT,pos1,pos2,k1,k2):
        
        #np.seterr(under="ignore")
        self.U=None            
        N = self.N
        self.dt=DT
        sigma=0.06
        print ("making gaussians in 4D")
        
        #self.psi4d = makeGaussian2P2D(N,L,pos1,pos2,k1,k2)        
        # make each particle psi with specific spin                                
        self.psi2dA = make2DGaussian(N,L, k1 , pos1,sigma)
        self.psi2dB = make2DGaussian(N,L, k2 , pos2,sigma)                                
        self.psi2dto4d()                                                
        
         # split step psi is calculated at each step DT step
        if self.V is None:                         
            self.V = np.full([N,N,N,N],1e-40)
                
        Lm=L*5.29177210903E-11
        
        DTsecs=DT*2.4188843265857E-17
        print("L meters ",Lm)
        print("DTsecs ",DTsecs)
        self.U = SplitStepMethod(self.V, (Lm, Lm,Lm, Lm), DTsecs) # units={'c': 1.0})
##
    def psi2dto4d(self):
        self.sumA=np.sum(self.psi2dA)
        self.sumB=np.sum(self.psi2dB)
        self.psi4d=np.einsum('jq,ki->ijkq', self.psi2dA, self.psi2dB)
        self.psi4d =(self.psi4d- np.transpose(self.psi4d, (1, 0, 3, 2)))/np.sqrt(2.0)

    def psi4dto2d(self):
        
        self.psi2dA = np.einsum('jqik->ij', self.psi4d)/self.sumB
        self.psi2dB= np.einsum('kiqj->ij', self.psi4d)/self.sumA
        
    def dostep(self):        
        self.psi4d = self.U(self.psi4d)
        
    def getProb(self):
        return np.einsum('ijkl->jk', np.abs(self.psi4d))
##        self.psi4dto2d()
##        p1 = np.abs(self.psi2dA)**2               
##        p2 = np.abs(self.psi2dB)**2        
##        return p1,p2
        
