# needs magpylib (2.3 pip3 )
# and easygui (but just for the choice box so u can easy remove it if dont want to install it)
import magpylib as magpy
from magpylib.source.magnet import Box
import numpy as np
from tqdm import tqdm
from plotmag import displaymagpsystem
from plotmag import plot_B_3D_arrows,animmagpysystem,displaymagpsystem
from easygui import choicebox
from threading import Thread

#fixed magnet parameters
#s = s /np.linalg.norm(s)
L = 2.  # space length in atomic units
N = 64
msg = "Select mode save file or load file and plot"
title = "B magnetic field"
choices = ["create B file npy Z  UP positive down negative",
           "create B file npy Z down positive top negative",
           "load B file npy and plot"]
choice = choicebox(msg, title, choices)

if choice==None:
    exit(0)
elif choice==choices[0]:
    Bdir=-1 # B positive top
elif choice==choices[1]:
    Bdir=1 # B positive down
elif choice==choices[2]:
    import tkinter.filedialog    
    filename = tkinter.filedialog.askopenfilename(filetypes=(        
        ("Archivos npy", "*.npy"),        
        ("Todos los archivos", "*.*") ))
    print(filename)
    if filename:            
        B=np.load(filename)
        print("B min",np.amin(B)," max ",np.amax(B))
        plot_B_3D_arrows(B, L,N)    
    exit(0)

def makeB3d(B,sg,L,N):    
    DX = (L/N)    
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
                B[x][y][z]=sg.getB([xx,yy,zz])                
                if zz>maxzz: maxzz=zz
                if zz<minzz: minzz=zz
                        
    print("minxx",minxx)
    print("minyy",minyy)
    print("minzz",minzz)
    print("maxxx",maxxx)
    print("maxyy",maxyy)
    print("maxzz",maxzz)
    print("Lmmhalf",Lmmhalf)
    

                        

    
borlenmm = 5.29177210903E-8 ## in milimeters
Lmm= L * borlenmm # space length in mm
B0=1500. # militeslas

### upper
mup = magpy.source.magnet.Box(mag=[0, -B0*Bdir, -B0*Bdir], dim=[Lmm*6, Lmm*2, Lmm*2],
                             pos=[0, 0, Lmm*2])
mup.rotate(45, [1,0,0])

# lower
mdown1 = magpy.source.magnet.Box(mag=[0, B0*Bdir,B0*Bdir], dim=[Lmm*6, Lmm*2, Lmm*2],
                             pos=[0, Lmm*1.4, -Lmm*2])
mdown1.rotate(45, [1,0,0])
mdown2 = magpy.source.magnet.Box(mag=[0, B0*Bdir,B0*Bdir], dim=[Lmm*6, Lmm*2, Lmm*2],
                             pos=[0, -Lmm*1.4, -Lmm*2])

mdown2.rotate(45, [1,0,0])
if Bdir== -1:
    mdown1.move([0,0,Lmm*4])
    mdown2.move([0,0,Lmm*4])
    mup.move([0,0,-Lmm*4])


c = magpy.Collection(mup,mdown1,mdown2)
displaymagpsystem(c,Lmm,10)
print("making 3d B field (will take  5 minutes)...")
B = np.zeros((N,N,N,3))
thmakeB = Thread(target=makeB3d, args=(B,c,L,N,))
thmakeB.start()
thmakeB.join()
if Bdir== 1:
    np.save("BOldPositiveDown.npy",B)
else:
    np.save("BOldPositiveUp.npy",B)
plot_B_3D_arrows(B, L,N)
print("end")            


