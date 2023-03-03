import numpy as np
#from dirac_splitstep import DiracSplitStepMethod
from auxfunctions import make2DGaussian,pauli4x4Matrixs,getEnergyEigenSpinors,initspinors3D
from splitstep import SplitStepMethod

class Dirac4D():
    
    def __init__(self,N):                
        self.M = 1.0 # Mass of the particle
        self.C =  137.036 # Speed of light
        self.HBAR = 1.0
        self.N=N
        self.dt=0.001
        self.U = None
        self.U_DAGGER = None
        self.E = None        
        self.psi4d=None
        self.psi2d=None
        self.V = None
        self.sumA=0
        self.sumB=0
        self.sumpsi2d=np.zeros([2,4],dtype=np.complex128)
        self.pauli4x, self.pauli4y, self.pauli4z = pauli4x4Matrixs()        
    def initialize(self,L,DT,pos1,pos2,k1,k2,spin1,spin2):
        
        #np.seterr(under="ignore")
        self.U=None            
        N = self.N
        self.dt=DT
        sigma=0.06
        print ("making gaussians in 4D")
        ones = np.ones([N,N], dtype=np.complex128)
        #self.psi4d = makeGaussian2P2D(N,L,pos1,pos2,k1,k2)
        
        # make each particle psi with specific spin        
                
        self.psi2d= np.zeros([2,4,N,N], dtype=np.complex128)
        
        kkk=np.array([k1,k2])
        
        pos=np.array([pos1,pos2])
        initial_spin=np.array([spin1,spin2])
        for i in range(2): 
            pos_eig1,pos_eig2, neg_eig1,neg_eig2 = getEnergyEigenSpinors(N,L,kkk[i],m=1.)            
            init_spinor  = initspinors3D(pos_eig1,pos_eig2,neg_eig1,neg_eig2, spin=initial_spin[i], bPositive=True)
            
            wavefunc = make2DGaussian(N,L, kkk[i] , pos[i],sigma)        
            self.psi2d[i] = wavefunc * np.multiply.outer(init_spinor, ones)
            
            #self.psi2d[i] /= np.linalg.norm(self.psi2d[i])
        
        
        
        #
        # get the 4d psi from the two 2d psi's
        
        self.psi4d= np.zeros([4,N,N,N,N], dtype=np.complex128)
        self.psi2dto4d()                        
        # test
##        print(self.getSpinExpecValue())
##        self.psi4dto2d()
##        print(self.getSpinExpecValue())
##        self.psi2dto4d()
##        self.psi4dto2d()
##        print(self.getSpinExpecValue())
##        exit(0)
        # end test
        
         # split step psi is calculated at each step DT step
        if self.V is None:                         
            self.V = np.full([N,N,N,N],1e-30)
        
        #self.U = DiracSplitStepMethod(self.V, (L,L,L,L), DT*2) # units={'c': 1.0})
        Lm=L*5.29177210903E-11
        DTm=DT*2.4188843265857E-17
        self.U = SplitStepMethod(self.V, (Lm, Lm,Lm, Lm), DTm) # units={'c': 1.0})
##
    def psi2dto4d(self):
                
        for i in range(4):
            self.sumpsi2d[0][i]=np.sum(self.psi2d[0][i])
            self.sumpsi2d[1][i]=np.sum(self.psi2d[1][i])
            
            self.psi4d[i]=np.einsum('jq,ki->ijkq', self.psi2d[0][i], self.psi2d[1][i])

    def psi4dto2d(self):
        for i in range(4):
            self.psi2d[0][i] = np.einsum('jqik->ij', self.psi4d[i]) /self.sumpsi2d[1][i]         
            self.psi2d[1][i]= np.einsum('kiqj->ij', self.psi4d[i]) / self.sumpsi2d[0][i]
        
    def dostep(self):
        for i in range(4):
            self.psi4d[i] = self.U(self.psi4d[i])
        
    def getProb(self):
        self.psi4dto2d()
            
        prob1 = sum([np.abs(self.psi2d[0][i])**2 for i in range(4)])
        prob2 = sum([np.abs(self.psi2d[1][i])**2 for i in range(4)])
        
        return prob1,prob2
    def getSpinExpecValue(self):
        sex=np.zeros([2,3],dtype=np.complex128)
        
        for i in range(2):
            sa = np.sum(np.abs(self.psi2d[i])**2)
            cp=np.conj(self.psi2d[i])
            sex[i][0] = np.sum(cp*np.einsum('ij,j...->i...', self.pauli4x, self.psi2d[i]))/sa
            sex[i][1] = np.sum(cp*np.einsum('ij,j...->i...', self.pauli4y, self.psi2d[i]))/sa
            sex[i][2] = np.sum(cp*np.einsum('ij,j...->i...', self.pauli4z, self.psi2d[i]))/sa
            sex[i] /= np.linalg.norm(sex[i])
                
        return sex
        
