import numpy as np
import matplotlib.pyplot as plt
import magpylib as magpy
from itertools import product, combinations
from mpl_toolkits.mplot3d import axes3d

from matplotlib import animation
def reduceNumberOfArrows(axes,B):
    Nx = len(axes[0])
    Ny = len(axes[1])
    Nz = len(axes[2])
    Barrows = np.zeros((Nx,Ny,Nz,3))
    for i in range (Nx):
        if i % 8 ==0:
            for j in range(Ny):
                if j % 8 ==0:
                    for k in range(Nz):
                        if k % 8 ==0:
                            Barrows[i][j][k]= B[i][j][k]
                            
    return Barrows



def plot_B_3D_arrows(B, L,N):    
    dim=3
    dx = L/N
    axes = [np.linspace(0, L, N)]*dim
    grid = np.meshgrid(*axes)
    
    maxB=np.amax(B)
    minB=np.amin(B)
    fig = plt.figure(figsize=(8,8),dpi = 100)

    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")
    ax.set_zlabel("z-axis")
    ax.set_title("Magnetic field ")

    ax.set_xlim(axes[0].min(), axes[0].max())
    ax.set_ylim(axes[1].min(), axes[1].max())
    ax.set_zlim(axes[2].min(), axes[2].max())
    ax.view_init(2,20)
    #ax.quiver(grid[0], grid[1], grid[2], J[:,:,:,1], J[:,:,:,0], J[:,:,:,2],length=0.1,color='darkblue')
    arrowlen=0.8/(maxB-minB)
    print("arrowlen",arrowlen)
    Br = reduceNumberOfArrows(axes,B)            
    ax.quiver(grid[0], grid[1], grid[2], Br[:,:,:,1], Br[:,:,:,0], Br[:,:,:,2],length=arrowlen,color='darkred')
    plt.tight_layout()
    def animation_func(i):            
        ax.view_init(1,i) #top or front
    anim = animation.FuncAnimation(fig, animation_func,                                   
                                   frames=360, interval=20, blit=False)
    plt.show()
def plotCube(ax, r):
    
    for s, e in combinations(np.array(list(product(r, r, r))), 2):
       if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s, e), color="green")
      
def displaymagpsystem(sg,Lmm,azimuth=30.):
    #define figure
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1,1,1, projection='3d')        
    ax.view_init(20,azimuth)
    plotCube(ax,[-Lmm/2., Lmm/2.])    
    figsg = magpy.displaySystem(sg,subplotAx=ax, direc=True)    
    ax.view_init(0,azimuth)
def animmagpysystem(sg,Lmm,topangle=0):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1,1,1, projection='3d')
    azimuth=0.
    ax.view_init(topangle,azimuth)
    plt.ion()    
    plotCube(ax,[-Lmm/2., Lmm/2.])
        
    
    while(azimuth<360.):
        
        ax.view_init(topangle,azimuth)
        
        magpy.displaySystem(sg,subplotAx=ax, direc=True)
        if azimuth==0:
            plt.pause(2)
        else: 
            plt.pause(0.2)
        azimuth+=2.
    plt.close()
