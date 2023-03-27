import numpy as np
import scipy.constants as const
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from numba import jit
from matplotlib import pyplot as plt


@jit(nopython=True)    
def applyrotmatrix4D( N: int , psi: np.ndarray, r4 : np.ndarray):
    # this rotations simulates a magnetic field
    sv= np.zeros(4,dtype=np.complex128)
    for x in range (N):
        for y in range (N):
            for a in range(N):
                for b in range(N):
                    sv[0]= psi[0][x][y][a][b]
                    sv[1]= psi[1][x][y][a][b]
                    sv[2]= psi[2][x][y][a][b]
                    sv[3]= psi[3][x][y][a][b]                                     
                    sv = np.dot(sv,r4[x][y][a][b])
                    psi[0][x][y][a][b]=sv[0]
                    psi[1][x][y][a][b]=sv[1]
                    psi[2][x][y][a][b]=sv[2]
                    psi[3][x][y][a][b]=sv[3]                               
                    
    

def arrayToRGBcmap(colors):
    npoints= len(colors)
    colormap = cm.jet
    normalize = mcolors.Normalize(vmin=np.min(colors), vmax=np.max(colors))
    s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
    rgb = s_map.to_rgba(colors).astype(np.float32)    
    return np.delete(rgb, 3,axis=1)  # reome alpha
        
def display4dpot(v4,i,umbral):
    array3d = np.abs(v4[:, :, :, i])    
    plt.rcParams["figure.figsize"] = [7.00, 3.50]
    plt.rcParams["figure.autolayout"] = True
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")
    ax.set_zlabel("z-axis")
    n =  len(v4[0][0][0])
    ax.set_xlim(0, n)
    ax.set_ylim(0,  n)
    ax.set_zlim(0, n)
    points =np.where( array3d > umbral)
    colors=arrayToRGBcmap(array3d[points])
    ax.scatter(points[0], points[1], points[2], c=colors, alpha=1)
    #x, y, z = np.where( array3d > umbral)
    #ax.scatter(x, y, z, c=z, alpha=1)
    
    plt.show()

def display3dpot(v3,umbral):    
    from matplotlib import pyplot as plt
    plt.rcParams["figure.figsize"] = [7.00, 3.50]
    plt.rcParams["figure.autolayout"] = True
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")
    ax.set_zlabel("z-axis")
    n =  len(v3[0][0])
    ax.set_xlim(0, n)
    ax.set_ylim(0,  n)
    ax.set_zlim(0, n)
    
    z, x, y = np.where( v3 > umbral)
    ax.scatter(x, y, z, c=z, alpha=1)
    plt.show()


@jit(nopython=True)
def pot2dtopot4d(n: int ,v2 : np.ndarray ,v4 : np.ndarray):    
    for x in range(n):            
        for y in range(n):
            for a in range(n):
                for b in range(n):                
                    v4[x][y][a][b]=v2[x][y]    
    
def makebarrier4d(n):
        v2=np.zeros([n,n])
        
        for x in range(n):            
            for y in range(12,16):                
                v2[x][y]=1e30
        v4=np.zeros([n,n,n,n])
        pot2dtopot4d(n,v2,v4)        
        return v4


def makeB2D(N,Bmaxmag,dirfield="x"):    

    B = np.zeros([3,N,N])                    
    vmag=-Bmaxmag            
    incmag=2.0 * Bmaxmag/self.N
        
    if dirfield == "x":            
        vdir=np.array([1.,0.])            
        for x in range (self.N):                
            vmag+=incmag
            vv = vdir*vmag
            for y in range(self.N):                                    
                B[0][x][y]= vv[0]
                B[1][x][y]= vv[1]
    
    elif dirfield == "y":                        
        vdir=np.array([0.,1.])                            
        for y in range (self.N):
            vmag+=incmag
            vv = vdir*vmag
            for x in range(self.N):                                    
                B[0][x][y]= vv[0]
                B[1][x][y]= vv[1]
    return B

def pauli2x2Matrixs():    
    σx = np.array([[0.0, 1.0], [1.0, 0.0]],np.complex128)
    σy = np.array([[0,-1j],[1j,0]],np.complex128)
    σz = np.array([[1,0],[0,-1]],np.complex128)
    return σx,σy,σz
def pauli4x4Matrixs():    
    p2x,p2y,p2z =  pauli2x2Matrixs()
    
    p4x = np.block([[p2x, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2x]])
    p4y = np.block([[p2y, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2y]])
    p4z = np.block([[p2z, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2z]])
    
    return p4x,p4y,p4z

def makeBpotential2D2(N,dt):
    # aplyed with
    #wave_func = np.einsum('ij...,j...->i...', exp_magnetic, wave_func).
    
    # load the B vectors stern gerlach 3D field in 3d and order axes for 2d
    # also set to 0  half of the field so it affects only to 1 particle
    #B = np.load("BgerlachZ32_2_5000.npy") *1.2

    
    #B = makeB2D(N,5000.,"y")
    B = np.load("BgerlachZ32_2_10k.npy") 
    B = B.transpose((3,1,2,0))        
    for a in range(N//2):
        for b in range(N):
            for c in range(N):
                    B[0][a][b][c]=0
                    B[1][a][b][c]=0
                    B[2][a][b][c]=0
    

    # exp(i q σ• B(r) Δt)
    # from https://en.wikipedia.org/wiki/Pauli_matrices#Exponential_of_a_Pauli_vector
    
    q = 1. #1.602176634e-19 #Coulombs
    a = q*np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)*dt
    d=np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    d[d == 0] = 1e-70    
    
    I = np.array([[1,0],[0,1]],np.complex128)        
    σx,σy,σz=pauli2x2Matrixs()
    ones=np.ones([N,N,N])
    I= np.multiply.outer(I, ones)
    σx= np.multiply.outer(σx, ones)
    σy= np.multiply.outer(σy, ones)
    σz= np.multiply.outer(σz, ones)
    
    ex = np.cos(a)*I - 1j* np.sin(a)* σx *B[0]/d
    ey = np.cos(a)*I + 1j* np.sin(a)* σy *B[1]/d
    ez = np.cos(a)*I + 1j* np.sin(a)* σz *B[2]/d
    pexp = ex+ey+ez

    pexp=pexp[:, :, :, :, 0]
        
    pexp[0][0] = pexp.real[0][0]/3. + (1j)*pexp.imag[0][0]
    pexp[1][1] = pexp.real[1][1]/3. + (1j)*pexp.imag[1][1]
        
   
    return pexp

def makeBpotential2D(N,dt):
    # aplyed with
    #wave_func = np.einsum('ij...,j...->i...', exp_magnetic, wave_func).
    
    # load the B vectors stern gerlach 3D field in 3d and order axes for 2d
    # also set to 0  half of the field so it affects only to 1 particle
    #B = np.load("BgerlachZ32_2_5000.npy") *1.2

    
    B = np.load("BgerlachZ32_2_10k.npy") 
    B = B.transpose((3,1,2,0))        
    for a in range(N//2):
        for b in range(N):
            for c in range(N):
                    B[0][a][b][c]=0
                    B[1][a][b][c]=0
                    B[2][a][b][c]=0
    

    # exp(i q σ• B(r) Δt)
    # from https://en.wikipedia.org/wiki/Pauli_matrices#Exponential_of_a_Pauli_vector
    
    q = 1. #1.602176634e-19 #Coulombs
    a = q*np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)*dt
    d=np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    d[d == 0] = 1e-70    
    
    nx, ny, nz = B/d
 
    a1=np.cos(a) + 1.0j*nz*np.sin(a)
    a2=1.0j*(nx - 1.0j*ny)*np.sin(a)
    b1=1.0j*(nx + 1.0j*ny)*np.sin(a)
    b2=np.cos(a) - 1.0j*nz*np.sin(a)
    
    a1=a1[:, :, 0]
    a2=a2[:, :, 0]
    b1=b1[:, :, 0]
    b2=b2[:, :, 0]
    
    
    pexp = np.array([[a1,a2],[b1,b2]])
##    plt.imshow(np.abs(pexp[1][1]))
##    plt.show()
##    exit(0)
    return pexp
       
def makeBpotential4D(N,dt):        
    B = np.load("BgerlachZ32_2_10k.npy") 
   
    B = B.transpose((3,1,2,0))        
##    for a in range(N//2):
##        for b in range(N):
##            for c in range(N):
##                    B[0][a][b][c]=0
##                    B[1][a][b][c]=0
##                    B[2][a][b][c]=0
        
    q = 1. #1.602176634e-19 #Coulombs
    # exp(i q σ• B(r) Δt)
    a = q*np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)*dt
    d=np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    d[d == 0] = 1e-70    
    
    I = np.array([[1,0],[0,1]],np.complex128)
    I = np.block([[I, np.zeros([2, 2])],
                [np.zeros([2, 2]), I]]) # 4x4 identity matrix
    σx,σy,σz=pauli4x4Matrixs()
    ones=np.ones([N,N,N])
    I= np.multiply.outer(I, ones)
    σx= np.multiply.outer(σx, ones)
    σy= np.multiply.outer(σy, ones)
    σz= np.multiply.outer(σz, ones)
    
    ex = np.cos(a)*I + 1j* np.sin(a)* σx *B[0]/d
    ey = np.cos(a)*I + 1j* np.sin(a)* σy *B[1]/d
    ez = np.cos(a)*I + 1j* np.sin(a)* σz *B[2]/d
    pexp = ex+ey+ez

    
    pexp=pexp[:, :, :, :, 0]
        
    pexp[0][0] = pexp.real[0][0]/3. + (1j)*pexp.imag[0][0]
    pexp[1][1] = pexp.real[1][1]/3. + (1j)*pexp.imag[1][1]
    print(np.shape(pexp))
    
        
    pexp4 = np.zeros([4,4,N,N,N,N],np.complex128)
    
    for i in range(4):
        for k in range(4):
            print("converting array 4x4:", i,k)
            p4=  np.zeros([N,N,N,N],np.complex128)
            pot2dtopot4d( N , pexp[i][k], p4)            
            pexp4[i][k]= np.copy(p4)
            
           
    print("potexp shape",np.shape(pexp4))
    
    return pexp4
        
##N=32    
##L=2.
##V = makeBpotential2D(N,0.001)
#i=0
#k=3
#umbral=np.amax(V[i][k])/2.
#display4dpot(V[i][k],N//2-1,umbral)

