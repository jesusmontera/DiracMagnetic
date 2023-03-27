# this class implements Pauli's with SplitStepMethod class
# and in my test's gives the same result as myPauli3D.py with eigen vectors
import numpy as np
from auxfunctions import make3DGaussian,spinDotB,pauli2x2Matrixs,makeBpotential,getBlochVector
from splitstep import SplitStepMethod
class myPauli3D():
    
    def __init__(self,N):                
        self.M = 1.0 # Mass of the particle
        self.C =  137.036 # Speed of light
        self.HBAR = 1.0
        self.N=N
        self.dt=0.001
        self.U = None        
        self.psi=None
        self.prob=[]
        self.spins=[]
        self.numframes=0                        
        self.pauli2x, self.pauli2y, self.pauli2z = pauli2x2Matrixs()
        self.initialspin=None
        self.exp_magnetic=None # B magnetic potential
        self.V = None
        self.Bmode ="pot"
    def clear(self):
        self.prob=[]
        self.psi=None
        self.spins=[]
        self.exp_magnetic=None
        self.numframes=0
            
        
    def getSpinExpecValue(self,psi):
        sa = np.sum(np.abs(psi)**2)
        cp=np.conj(psi)
        sx = np.sum(cp*np.einsum('ij,j...->i...', self.pauli2x, psi))/sa
        sy = np.sum(cp*np.einsum('ij,j...->i...', self.pauli2y, psi))/sa
        sz = np.sum(cp*np.einsum('ij,j...->i...', self.pauli2z, psi))/sa        
        sexpec=np.real(np.array([sx,sy,sz],np.complex128))
        n = np.linalg.norm(sexpec)
        if n!=0:
            sexpec/=n            
        return sexpec
    
    def get_energies(self,momenta):
        C=self.C
        px, py, pz = momenta
        omega = np.sqrt(self.M**2*C**2 + px**2 + py**2 + pz**2)
        return np.array([C*omega, C*omega])

    
                        
    def initAnimation(self,L,DT,k0, pos0,initial_spin=None, B=None,Bmode="pot"):
        
        np.seterr(under="ignore")
        self.prob= []
        self.spins=[]
        self.numframes=0
        self.initialspin = initial_spin
        N = self.N
        self.dt=DT
        self.Bmode=Bmode
        if B is not None and Bmode=="pot":                        
            self.exp_magnetic = makeBpotential( DT, B.transpose((3,0,1,2)) )
                    
        ones = np.ones([N,N,N], dtype=np.complex128)
        
        sigma=0.06
        wavefunc = make3DGaussian(N,L, k0 , pos0,sigma)
        #print("wf sum",np.sum(np.abs(wavefunc)**2))        
        self.psi = wavefunc * np.multiply.outer(initial_spin, ones)
        n = np.linalg.norm(self.psi)
        if n != 0:
            self.psi /= n
        
        
        Lm=L*5.29177210903E-11
        DTm=DT*2.4188843265857E-17
        if self.V is None:            
            self.V= np.zeros([N,N,N])
            
        self.U = SplitStepMethod(self.V, (Lm, Lm, Lm), DTm) # units={'c': 1.0})

            
        
        
        
        
        
    def saveFileProb (self, file):
        np.save(file, self.prob,allow_pickle=True)
        spinfile = file.replace(".npy","spin.npy")
        np.save(spinfile, self.spins,allow_pickle=True)
        
    def loadFileProb (self, file):
        self.prob = np.load(file,allow_pickle=True)
        spinfile = file.replace(".npy","spin.npy")
        self.spins = np.load(spinfile,allow_pickle=True)
        self.numframes=len(self.prob)
        return self.numframes
    def isLoaded(self):
        if self.prob==[]:
            return False
        else: return True
    def getProbability(self,i):        
        return self.prob[i]
    def getSpinbloch(self,i):

        posspin=self.spins[i*2]
        spinbloch=self.spins[i*2+1]                            
        return posspin, spinbloch


    def doAnimFrame(self,Bmagnetic=None,Usteps=1):                

        if self.numframes>0:
            for _ in range(Usteps):
                for i in range(2):
                    self.psi[i] = self.U(self.psi[i])                

                if Bmagnetic is not None:
                    
                    sbloch= self.getSpinExpecValue(self.psi)
                    if self.Bmode=="pot": #with B as potential, spin rotates so looks ok
                        self.psi = np.einsum('ij...,j...->i...', self.exp_magnetic, self.psi)
                    else:  # with B as dot, spin don't rotate so is a botch
                        spinDotB(self.N, self.dt, self.psi,Bmagnetic, sbloch)
                n =np.linalg.norm(self.psi)                
                if n!=0:
                    self.psi /= n
        
        # save prob in array
        self.prob.append( sum([np.abs(self.psi[i])**2 for i in range(2)]))

        # and save spin  in array (it's position is at  max prob)
        p = self.prob[self.numframes]        
        posmax = np.unravel_index(p.argmax(), p.shape) # max index  from prob is spin pos
        self.spins.append(np.array(posmax))                
        if Bmagnetic is None or self.numframes==0:
            sbloch= self.getSpinExpecValue(self.psi)        
        self.spins.append(sbloch)
       # print("frame ",self.numframes, " spin",np.round(sbloch,3))    
        self.numframes+=1

    
        

       
       
    
