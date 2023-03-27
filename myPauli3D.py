import numpy as np
from auxfunctions import make3DGaussian,spinDotB,pauli2x2Matrixs,makeBpotential,getBlochVector,applyDotMatrix3D
from splitstep import SplitStepMethod
class myPauli3D():
    
    def __init__(self,N):                
        self.M = 1.0 # Mass of the particle
        self.C =  137.036 # Speed of light
        self.HBAR = 1.0
        self.N=N
        self.dt=0.001
        self.U = None
        self.U_DAGGER = None
        self.E = None        
        self.psi=None
        self.prob=[]
        self.spins=[]
        self.numframes=0                        
        self.pauli2x, self.pauli2y, self.pauli2z = pauli2x2Matrixs()
        self.initialspin=None
        self.exp_magnetic=None # B magnetic potential
        self.Bmode ="pot"
        self.V = None
        self.bUseSplit=False # set manually because is not implement in graphic interface
    def clear(self):
        self.prob=[]
        self.U = None
        self.U_DAGGER = None
        self.E = None        
        self.psi=None
        self.V=None
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
##        n = np.linalg.norm(sexpec)
##        if n!=0:
##            sexpec/=n            
        return sexpec
    
    def get_energies(self,momenta):
        C=self.C
        px, py, pz = momenta
        omega = np.sqrt(self.M**2*C**2 + px**2 + py**2 + pz**2)
        return np.array([C*omega, C*omega])
    def get_eigenvectors(self,L,k):
        """
        For the given momenta find the corresponding spinor
        energy eigenvectors that solves the time-independent
        free-particle Dirac equation. This returns a 2x2 matrix
        where each row is the eigenvector.
        """
        px, py, pz = 2.0*k[0]*np.pi/L, 2.0*k[1]*np.pi/L, 2.0*k[2]*np.pi/L
        #px, py, pz = momenta
        p2 = px**2 + py**2 + pz**2
        p = np.sqrt(p2)
        mc = self.M*self.C        
        omega = np.sqrt(mc*mc + p2) # Corresponds to E/c

        # Temporary variable used for
        # some denominators in the negative energy solutions
        den1 = p*np.sqrt((mc - omega)**2 + p2)
        # Used for denominators in the positive energy solutions
        den2 = p*np.sqrt((mc + omega)**2 + p2)

        
        # Positive energy solutions
        
        pos_eig1 = [pz*(mc + omega)/den2, 
                    (mc*px - 1.0j*mc*py + (px - 1.0j*py)*omega)/den2]    
        pos_eig2 = [(mc*px + 1.0j*mc*py + (px + 1.0j*py)*omega)/den2,
                    -pz*(mc + omega)/den2]

        
        return np.array([pos_eig1, pos_eig2])

    @staticmethod
    def initspinors3D(pos_eig1,pos_eig2, spin=None):    
        
        if spin is None:
            # if no initial spin is requiered  it will point to Z up
            spinorsketpos = np.array([pos_eig1[1], pos_eig1[0]])            
            print("init spinors from momentum with Z up")
            
        else:        
            # if initial spin is requiered  make spinors for that specific spin
            # convert spinors to ket form to merge spin with momentum
            print("init spinors from specific spin" , np.round(spin,3),
                  "\n\tinitial spin Bloch = ",np.round(getBlochVector(spin),3))
            c1 = pos_eig1 @ np.array([spin[0], spin[1]]) # inner product
            c2 = pos_eig2 @ np.array([spin[0], spin[1]])
            spinorsketpos = c1 * np.conj(pos_eig1) + c2 * np.conj(pos_eig2)
            

        n = np.linalg.norm(spinorsketpos)        
        if n!=0:
            spinorsketpos = spinorsketpos / n
        
        init_spinor = [spinorsketpos[0],spinorsketpos[1]]
        
        return init_spinor
    @staticmethod
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
        # Positive energy eigenstates in bra form
                
        pos_eig1 = [pz*(mc + omega)/den2, 
                    (mc*px - 1.0j*mc*py + (px - 1.0j*py)*omega)/den2]
        pos_eig2 = [(mc*px + 1.0j*mc*py + (px + 1.0j*py)*omega)/den2,
                    -pz*(mc + omega)/den2]
            
        return pos_eig1,pos_eig2
                        
    def initAnimation(self,L,DT,k0, pos0,initial_spin=None, B=None,Bmode="pot"):
        
        np.seterr(under="ignore")
        self.clear()
        self.initialspin = initial_spin
        N = self.N
        self.dt=DT
        self.Bmode=Bmode
        if B is not None:# and Bmode=="pot":
            #B transposed is [3,N,N,N]
            self.exp_magnetic = makeBpotential( DT,  B.transpose((3,0,1,2)) )
            
                    
        ones = np.ones([N,N,N], dtype=np.complex128)
        pos_eig1,pos_eig2 = myPauli3D.getEnergyEigenSpinors(N,L,k0,m=1.)        
        init_spinor  = myPauli3D.initspinors3D(pos_eig1,pos_eig2, spin=initial_spin)        
        sigma=0.06
        wavefunc = make3DGaussian(N,L, k0 , pos0,sigma)
        #print("wf sum",np.sum(np.abs(wavefunc)**2))        

        self.psi = wavefunc * np.multiply.outer(init_spinor, ones)
        #or replace with above with(but inital spin is required )
        #self.psi = wavefunc * np.multiply.outer(initial_spin, ones)
        
        n = np.linalg.norm(self.psi)
        if n != 0:
            self.psi /= n
        
        
        # free periodic (psi can be calculate at any time, but because
        # there can be a magnetic field(B) is calculated in DT time steps

        
        X, Y, Z = np.meshgrid(L*np.linspace(-0.5, 0.5 - 1.0/N, N),
                              L*np.linspace(-0.5, 0.5 - 1.0/N, N),
                              L*np.linspace(-0.5, 0.5 - 1.0/N, N))

        if self.bUseSplit:
            Lm=L*5.29177210903E-11
            DTm=DT*2.4188843265857E-17
            self.V = np.zeros([N,N,N])        
            self.U = SplitStepMethod(self.V, (Lm, Lm,Lm), DTm) 
        else:
            f = np.fft.fftfreq(N)
            f[0] = 1e-60
            P = np.pi * f * N / L # Momenta in 1D
            PX, PY, PZ = np.meshgrid(P, P, P) # Momenta in the x y z directions
            
            self.E = self.get_energies([PX, PY, PZ])
    ##        self.U =self.get_eigenvectors(L,[PX, PY, PZ])
                                                                   
            self.U = np.array([np.multiply.outer(pos_eig1, ones),
                               np.multiply.outer(pos_eig2, ones)])

            ind = [i for i in range(3 + 2)]
            ind[0], ind[1] = ind[1], ind[0]
            self.U_DAGGER = np.conj(np.transpose(self.U, ind))
            
        
        
        
        
        
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

    def dostep(self,psi,dt):
        # go to momentun space
        psi = np.array([np.fft.fftn(psi[i]) for i in range(2)])
        #wavefunction in terms of spinor energy eigenvectors.
        psi = np.einsum('jk...,k...->j...', self.U, psi)
        #do time step
        psi = np.exp(-1.0j*dt*self.E/self.HBAR)*psi
        psi = np.einsum('ij...,j...->i...', self.U_DAGGER, psi)
        for i in range(2):
            psi[i] = np.fft.ifftn(psi[i])
        return psi

    def doAnimFrame(self,Bmagnetic=None,Usteps=1):                

        if self.numframes>0:
            for _ in range(Usteps):
                if self.bUseSplit:
                    for i in range(2):                    
                        self.psi[i] = self.U(self.psi[i])
                else:
                    self.psi = self.dostep(self.psi,self.dt *2. )
                
                if Bmagnetic is not None:
                    
                    sbloch= self.getSpinExpecValue(self.psi)
                    if self.Bmode=="pot": #with B as potential, spin rotates so looks ok
                        self.psi = np.einsum('ij...,j...->i...', self.exp_magnetic, self.psi)
                    else:  # with B as dot, spin don't rotate so is a botch
                        bv= np.zeros([2,2],dtype=np.complex128)
                        #print(np.shape(self.exp_magnetic))
                        applyDotMatrix3D( self.N , self.psi,self.exp_magnetic,bv)                                            
                        #spinDotB(self.N, self.dt, self.psi,Bmagnetic, sbloch)
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

    
        

       
       
    
