import magpylib as magpy
from magpylib.source.magnet import Box
import numpy as np
from tqdm import tqdm
from plotmag import displaymagpsystem
from plotmag import plot_B_3D_arrows,animmagpysystem,displaymagpsystem
from easygui import choicebox
#fixed magnet parameters
#s = s /np.linalg.norm(s)
L = 2.  # space length in atomic units
N = 64
msg = "Select mode save file or load file and plot"
title = "B magnetic field"
choices = ["create B file npy",        
           "load B file npy and plot"]
choice = choicebox(msg, title, choices)
if choice==choices[1]:
    B=np.load("Bgerlach.npy")
    print("B min",np.amin(B)," max ",np.amax(B))
    plot_B_3D_arrows(B, L,N)
    
    exit(0)

def makeB3d(sg,L,N):    
    DX = (L/N)
    B = np.zeros((N,N,N,3))
    Lmmhalf = Lmm/2.
    NtoLmm= DX * borlenmm
    maxxx,maxyy,maxzz=-1E10,-1E10,-1E10
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
borlenmm = 5.29177210903E-8 ## in milimeters
Lmm= L * borlenmm # space length in mm
B0=1000 # militeslas
### upper
mup = magpy.source.magnet.Box(mag=[0, -B0, -B0], dim=[Lmm*6, Lmm*2, Lmm*2],
                             pos=[0, 0, Lmm*2])
mup.rotate(45, [1,0,0])

# lower
mdown1 = magpy.source.magnet.Box(mag=[0, B0,B0], dim=[Lmm*6, Lmm*2, Lmm*2],
                             pos=[0, Lmm*1.4, -Lmm*2])
mdown1.rotate(45, [1,0,0])
mdown2 = magpy.source.magnet.Box(mag=[0, B0,B0], dim=[Lmm*6, Lmm*2, Lmm*2],
                             pos=[0, -Lmm*1.4, -Lmm*2])

mdown2.rotate(45, [1,0,0])



c = magpy.Collection(mup,mdown1,mdown2)
displaymagpsystem(c,Lmm)
print("making 3d B field...")
Bfield = makeB3d(c,L,N)
np.save("Bgerlach.npy",Bfield) 
plot_B_3D_arrows(Bfield, L,N)

print("end")                                    

