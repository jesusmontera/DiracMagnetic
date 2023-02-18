import numpy as np
from dirac_splitstep import DiracSplitStepMethod
from auxfunctions import getBlochVector,make3DGaussian,getPsiBlochSpin,spinDotBdirac

class myDirac3D():
    
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
        self.bSplit=False
        self.V=np.zeros((N,N,N))
    def clear(self):
        self.prob=[]
        self.spins=[]
        self.numframes=0
   
    
    def get_energies(self,momenta):
        C=self.C
        px, py, pz = momenta
        omega = np.sqrt(self.M**2*C**2 + px**2 + py**2 + pz**2)
        return np.array([-C*omega, -C*omega, C*omega, C*omega])
    def get_eigenvectors(self,k):
        """
        For the given momenta find the corresponding spinor
        energy eigenvectors that solves the time-independent
        free-particle Dirac equation. This returns a 4x4 matrix
        where each row is the eigenvector.
        """
        px, py, pz = 2.0*k[0]*np.pi/L, 2.0*k[1]*np.pi/L, 2.0*k[2]*np.pi/L
        #px, py, pz = momenta
        p2 = px**2 + py**2 + pz**2
        p = np.sqrt(p2)
        mc = self.M*self.C
        zeros = np.zeros(3*[self.N], dtype=np.complex128)
        omega = np.sqrt(mc*mc + p2) # Corresponds to E/c

        # Temporary variable used for
        # some denominators in the negative energy solutions
        den1 = p*np.sqrt((mc - omega)**2 + p2)
        # Used for denominators in the positive energy solutions
        den2 = p*np.sqrt((mc + omega)**2 + p2)


        # Negative energy solutions
          
        neg_eig1 = [pz*(mc - omega)/den1, 
                    (mc*px - 1.0j*mc*py - (px - 1.0j*py)*omega)/den1,
                    p2/den1, zeros]

        neg_eig2 = [(mc*px + 1.0j*mc*py + (-px - 1.0j*py)*omega)/den1,
                    -pz*(mc - omega)/den1, zeros, p2/den1]

        # Positive energy solutions
        
        pos_eig1 = [pz*(mc + omega)/den2, 
                    (mc*px - 1.0j*mc*py + (px - 1.0j*py)*omega)/den2,
                    p2/den2, zeros]    
        pos_eig2 = [(mc*px + 1.0j*mc*py + (px + 1.0j*py)*omega)/den2,
                    -pz*(mc + omega)/den2, zeros, p2/den2]
        
        
        return np.array([neg_eig1, neg_eig2, pos_eig1, pos_eig2])

    @staticmethod
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
        #spin0spinors = spinQubitfromspinor(spinorsketpos,bPositive)
        print(sinfo)            
        init_spinor = [spinorsketpos[0],spinorsketpos[1],spinorsketpos[2],spinorsketpos[3] ]
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
                        
    def initAnimation(self,bSplit,L,DT,k0, pos0,initial_spin=None):
        
        np.seterr(under="ignore")
        self.prob= []
        self.spins=[]
        self.numframes=0
        
        N = self.N
        self.dt=DT
        self.bSplit=bSplit
        ones = np.ones([N,N,N], dtype=np.complex128)
        pos_eig1,pos_eig2, neg_eig1,neg_eig2 = myDirac3D.getEnergyEigenSpinors(N,L,k0,m=1.)
        
        if bSplit:
            # split step psi is calculated at each step DT step
            self.U = DiracSplitStepMethod(self.V, (L, L, L), DT*2) # units={'c': 1.0})                        
            
        else:
            # free periodic (psi can be calculate at any time, but because
            # there can be a magnetic field(B) is calculated in DT time steps
            X, Y, Z = np.meshgrid(L*np.linspace(-0.5, 0.5 - 1.0/N, N),
                                  L*np.linspace(-0.5, 0.5 - 1.0/N, N),
                                  L*np.linspace(-0.5, 0.5 - 1.0/N, N))

            f = np.fft.fftfreq(N)
            f[0] = 1e-60
            P = np.pi * f * N / L # Momenta in 1D
            PX, PY, PZ = np.meshgrid(P, P, P) # Momenta in the x y z directions
            
            self.E = self.get_energies([PX, PY, PZ])
            #self.U = self.get_eigenvectors([PX, PY, PZ]) # 4x4 matrix containing the eigenvectors
            self.U = np.array([np.multiply.outer(neg_eig1, ones),
                               np.multiply.outer(neg_eig2, ones),
                               np.multiply.outer(pos_eig1, ones),
                               np.multiply.outer(pos_eig2, ones)])

            ind = [i for i in range(3 + 2)]
            ind[0], ind[1] = ind[1], ind[0]
            self.U_DAGGER = np.conj(np.transpose(self.U, ind))
            
        sigma=0.06
        wavefunc = make3DGaussian(N,L, k0 , pos0,sigma)
        initial_spin=None #np.array([0.70710678+0.j, 0.70710678+0.j],np.complex128)                        
        init_spinor  = myDirac3D.initspinors3D(pos_eig1,pos_eig2,neg_eig1,neg_eig2, spin=initial_spin, bPositive=True) 
        self.psi = wavefunc * np.multiply.outer(init_spinor, ones)
        self.psi /= np.linalg.norm(self.psi)        
        #print("wf sum",np.sum(np.abs(wavefunc)**2))        
        
        
        
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
        psi = np.array([np.fft.fftn(psi[i]) for i in range(4)])
        #wavefunction in terms of spinor energy eigenvectors.
        psi = np.einsum('jk...,k...->j...', self.U, psi)
        #do time step
        psi = np.exp(-1.0j*dt*self.E/self.HBAR)*psi
        psi = np.einsum('ij...,j...->i...', self.U_DAGGER, psi)
        for i in range(4):
            psi[i] = np.fft.ifftn(psi[i])
        return psi

    def doAnimFrame(self,Bmagnetic=None,modeSpin=0,Usteps=1):                

        
        for _ in range(Usteps):
            if self.numframes>0:
                if self.bSplit:
                    self.psi = self.U(self.psi)
                else:
                    self.psi = self.dostep(self.psi,self.dt *2. )

                if Bmagnetic is not None:                
                    spinpsibloch = getPsiBlochSpin(self.psi,modeSpin)
                    spinDotBdirac(self.N, self.dt, self.psi,Bmagnetic, spinpsibloch)
                self.psi /= np.linalg.norm(self.psi)
        
        # save prob in array
        self.prob.append( sum([np.abs(self.psi[i])**2 for i in range(4)]))

        # and save spin  in array (it's position is at  max prob)
        p = self.prob[self.numframes]
        posmax = np.unravel_index(p.argmax(), p.shape) # max index  from prob is spin pos
        self.spins.append(np.array(posmax))        
        spinpsibloch = getPsiBlochSpin(self.psi, modeSpin)            
        self.spins.append(spinpsibloch)        
        self.numframes+=1

    
        

       
       
    
