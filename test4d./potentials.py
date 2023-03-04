import numpy as np
import scipy.constants as const

def makeYBarrier(N,Lmeters):
    X = Lmeters*np.linspace(-0.5, 0.5 - 1.0/N, N)
    X1, X2, Y1, Y2 = np.meshgrid(X, X, X, X)
    R=np.abs(Y1) + 1e-60    #Y1
    V4d = const.e**2/(4.0*np.pi*const.epsilon_0*R)
    
    return V4d

        
        
    
