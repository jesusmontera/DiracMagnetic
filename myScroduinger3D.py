import numpy as np
from splitstep import SplitStepMethod
from auxfunctions import getBlochVector,make3DGaussian
from numba import jit

class mySchroduinger3D():
    def __init__(self,N):                
        self.N=N        
        #self.L_METERS =         
        self.U=None
        self.V = np.zeros((N,N,N))        
        self.data=None
        self.prob=[]
        self.blochspin=np.array([1.,0.,0.])
        self.spin=[]
        self.numframes=0
        
    def initAnimation(self,L,DT,k0, POS0, initial_spin=None):
        # p si momentum
        np.seterr(under="ignore")
        
        #print("schroduinger DT in seconds",DT*2.4188843265857E-17)
        #print("schroduinger L in meters",L*5.29177210903E-11)
        Lm=L*5.29177210903E-11
        DTm=DT*2.4188843265857E-17  
        
        self.U = SplitStepMethod(self.V, (Lm, Lm, Lm), DTm * 2.) # 4600 fix to adjust time       
        wavefunc = make3DGaussian(N=self.N,L=L, k=k0 , pos = POS0,sigma=0.07)
       
       
        self.data = {'psi': wavefunc}
        self.data['psi'] /= np.linalg.norm(self.data['psi'])
        print("wf sum",np.sum(np.abs(self.data['psi'] )**2))        
        self.prob= []
        self.blochspin= initial_spin       
        self.numframes=0
        
    def saveFileProb (self, file):
        np.save(file, self.prob,allow_pickle=True)
        spinfile = file.replace(".npy","spin.npy")
        np.save(spinfile, self.spin,allow_pickle=True)
        
    def loadFileProb (self, file):
        self.prob = np.load(file,allow_pickle=True)
        spinfile = file.replace(".npy","spin.npy")
        self.spin = np.load(spinfile,allow_pickle=True)
        self.numframes=len(self.prob)
        return self.numframes
    def clear(self):
        self.prob=[]
        self.spin=[]
        self.numframes=0
    def isLoaded(self):
        if self.prob==[]:
            return False
        else: return True
    def getProbability(self,i):        
        return self.prob[i]
    def getSpinbloch(self,i):
        posspin=self.spin[i*2]        
        spinbloch = self.spin[i*2+1]                
        return posspin, spinbloch

    @staticmethod
    @jit(nopython=True)
    def spinDotB(N: int, DT: float, wf: np.ndarray,B: np.ndarray, spinbloch: np.ndarray):
        maxdd=-1e30
        mindd=1e30
        for x in range(N):
            for y in range(N):
                for z in range(N):                    
                    R = wf[x][y][z].real
                    I = wf[x][y][z].imag                                        
                    # imag
                    dd = spinbloch.dot(B[x][y][z]) #* prob #(/maxprob)
                    if dd > maxdd:  maxdd=dd
                    if dd < mindd:  mindd=dd
                        
                    Iinc =  dd * R * DT
                    #real
                    Rdec = dd * I * DT 
                    inc= -Rdec + Iinc *1j
                    #apply B contribution
                    wf[x][y][z] += inc
        #print("spinDotB mindd = ",mindd,"maxdd",maxdd)
    def doAnimFrame(self,dt=0.01,Bmagnetic=None,Usteps=1):        
                    
        for _ in range(Usteps):        
            self.data['psi'] = self.U(self.data['psi'])
            if Bmagnetic is not None:                
                self.spinDotB(self.N, dt, self.data['psi'],Bmagnetic, self.blochspin)
            self.data['psi'] /= np.linalg.norm(self.data['psi'])
                
        # save prob in array
        self.prob.append(np.abs(self.data['psi']))
        # and save spin  in array (it's position is at  max prob)
        p = self.prob[self.numframes]
        posmax = np.unravel_index(p.argmax(), p.shape) # max index  from prob is spin pos
        self.spin.append(np.array(posmax))
        self.spin.append(self.blochspin)
        self.numframes+=1
        
    
        

       
       
    
