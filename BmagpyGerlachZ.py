# create a Stern-Gerlach B-field approximation with 3 cuboid magnets
from itertools import product, combinations
from magpylib.source.magnet import Box
from magpylib import Collection,displaySystem


import matplotlib.pyplot as plt
import numpy as np

def plotCube(ax, r):    
    for s, e in combinations(np.array(list(product(r, r, r))), 2):
       if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s, e), color="green")

def displaymagpsystem(sg,Lmm,azimuth=0.):
    #define figure
    fig = plt.figure(figsize=(6, 6))    
    
    ax = fig.add_subplot(1,1,1, projection='3d')        
    ax.view_init(0,azimuth)
    plotCube(ax,[-Lmm/2., Lmm/2.])        
    displaySystem(sg,subplotAx=ax, direc=True)
##    plt.pause(3)
##    plt.close()
    return fig
    

def MakeBmagpyZ(N=64,L_au=2.0,B0=2000.0, Bdir=-1):
    # B0 miliTeslas
    # Bdir 1 or -1
# define constants
    a0_mm = 5.29177210903e-8  ## Bohr length in millimeters
    L_mm = L_au * a0_mm  # space side length in millimeters
    B1 = B0 / (2 ** 0.5)

    # create magnetic system
    mag_dim = [L_mm * 6, L_mm * 2, L_mm * 2]  # cuboids dimensions
    mag_upper = Box(
        mag=[0, -B1, -B1],
        dim=mag_dim,
        pos=[0, 0, L_mm * 2],
    )
    mag_upper.rotate(45, [1, 0, 0])

    mag_lower1 = Box(
        mag=[0, B1, B1],
        dim=mag_dim,
        pos=[0, L_mm * 1.4, -L_mm * 2],
    )
    mag_lower1.rotate(45, [1, 0, 0])
    mag_lower2 = Box(
        mag=[0, B1, B1],
        dim=mag_dim,
        pos=[0, -L_mm * 1.4, -L_mm * 2],
    )
    mag_lower2.rotate(45, [1, 0, 0])

    coll = Collection(
        mag_upper, mag_lower1, mag_lower2)


    # flip system arount x-axis if field direction `Bdir` is negative
    if Bdir == -1:
        coll.rotate(180, [1, 0, 0], anchor=(0, 0, 0))

    

    # compute B_array
    ls = np.linspace(
        -L_mm / 2, L_mm / 2, N, endpoint=False
    )  # add endpoint=False if you want same behavior as in your example
    points = np.array([[x, y, z] for x in ls for y in ls for z in ls])
    print("making magpy B array in z ...wait a few seconds")
    
    B_array = coll.getB(points).reshape(N, N, N, 3)
    figsg = displaymagpsystem(coll,L_mm,10)
    
    
    return B_array
