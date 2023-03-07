import magpylib as magpy# version 2.3.0b0
import numpy as np

from tqdm import tqdm
from plotmag import plot_B_3D_arrows,displaymagpsystem
from easygui import choicebox

borlenmm = 5.29177210903E-8 ## in milimeters
L = 2.  # space length in atomic units
N = 64
Lmm= L * borlenmm # space length in mm
B0=1000 # militeslas

msg = "Select mode save file or load file and plot"
title = "B magnetic field"
choices = ["create B file npy",        
           "load B file npy and plot"]
choice = choicebox(msg, title, choices)
if choice==choices[1]:
    B=np.load("Bmio.npy")
    print("B min",np.amin(B)," max ",np.amax(B))
    plot_B_3D_arrows(B, L,N)
    
    exit(0)


magnet1 = magpy.source.magnet.Box(mag=[0, 0, B0], dim=[Lmm*2, Lmm*2, Lmm], pos=[0, 0, -Lmm*2])
magnet2 = magpy.source.magnet.Box(mag=[0, 0, -B0], dim=[Lmm*2, Lmm*2, Lmm], pos=[0, 0, Lmm*2])
# collect magnets
sg = magpy.Collection(magnet1,magnet2)

def makeB3d(sg,L,N):    
    DX = (L/N)
    B = np.zeros((N,N,N,3))
    Lmmhalf = Lmm/2.
    NtoLmm= DX * borlenmm
    maxxx,maxyy,maxzz=-1E10.,-1E10,-1E10.
    minxx,minyy,minzz=1E10,-1E10,-1E10
    for x in tqdm(range(N)):        
        xx= x * NtoLmm  -  Lmmhalf
        if xx>maxxx: maxxx=xx
        if xx<minxx: minxx=xx
        for y in range(N):
            yy= y * NtoLmm - Lmmhalf
            if yy>maxyy: maxyy=yy
            if yy<minyy: minyy=yy
            for z in range(N):
                zz= z * NtoLmm - Lmmhalf 
                B[x][y][z] = sg.getB([xx,yy,zz])
                if zz>maxzz: maxzz=zz
                if zz<minzz: minzz=zz
                        
    print("minxx",minxx)
    print("minyy",minyy)
    print("minzz",minzz)
    print("maxxx",maxxx)
    print("maxyy",maxyy)
    print("maxzz",maxzz)
    print("Lmmhalf",Lmmhalf)
    return B

#magpy.displaySystem(sg)
displaymagpsystem(sg,Lmm)
print("making 3d B field...")

Bfield = makeB3d(sg,L,N)
np.save("Bmio.npy",Bfield) 
plot_B_3D_arrows(Bfield, L,N)

print("end")                                    





