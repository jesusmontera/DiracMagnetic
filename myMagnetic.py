import numpy as np

class magneticlass():
    def __init__(self,N):
                
        self.N = N                    
        self.B = None
        self.FIELD_X=0
        self.FIELD_Y=1
        self.FIELD_Z=2
        self.dirfield=0 # x y o z
        self.bInvert=False 
        self.Bmaxmag=10.
        self.B =None # magnetuic filed N*N*N*3

        # list with two 3d vectors for draw main arrow
        # (first  position and second dir vector )
        self.mainArrow=[]
        self.btext=""
    def clear(self):
        self.btext=""
        self.mainArrow=[]
        self.B =None
    def load_npyFile(self, sfilepath):
        self.mainArrow=[]
        self.btext="magpy"
        self.B = np.load(sfilepath)
        

                    
            
        
    def _makeMainArrow(self):
        
        if self.bInvert:
            if self.dirfield == self.FIELD_X:
                vpos=np.array([self.N//2+12,0,0])
                vdir=np.array([-1.,0.,0.])
            elif self.dirfield == self.FIELD_Y:
                vpos=np.array([0,self.N//2+12,0])
                vdir=np.array([0,-1,0])
                
            elif self.dirfield == self.FIELD_Z:
                vpos=np.array([0,-self.N//2-6,0])
                vdir=np.array([0,0,-1])
        else:
            if self.dirfield == self.FIELD_X:
                vpos=np.array([-self.N//2-12,0,0])
                vdir=np.array([1,0,0])
            elif self.dirfield == self.FIELD_Y:
                vpos=np.array([0,-self.N//2-12,0])
                vdir=np.array([0,1,0])
            elif self.dirfield == self.FIELD_Z:
                vpos=np.array([0,self.N//2+6,0])
                vdir=np.array([0,0,1.])
        self.mainArrow.append(vpos)
        self.mainArrow.append(vdir)
                
            
            
    def makeB(self,dirfield,Bmaxmag,bInvert):
        self.btext=""
        self.B = np.zeros((self.N,self.N,self.N,3))
        vmag=0.
        self.dirfield=dirfield
        incmag=Bmaxmag/self.N
        self.Bmaxmag=Bmaxmag
        self.bInvert=bInvert
        if dirfield == self.FIELD_X:
            vdir=np.array([1.,0.,0.])
            for x in range (self.N):
                vmag+=incmag
                if bInvert:
                    vv = vdir*(Bmaxmag-vmag)
                else:
                    vv = vdir*vmag
                for y in range(self.N):                    
                    for z in range(self.N):                        
                        self.B[x][y][z]= vv
        
        elif dirfield == self.FIELD_Y:
            vdir=np.array([0.,1.,0.])
            for y in range (self.N):
                vmag+=incmag
                if bInvert:
                    vv = vdir*(Bmaxmag-vmag)
                else:
                    vv = vdir*vmag
                for x in range(self.N):                    
                    for z in range(self.N):                        
                        self.B[x][y][z]= vv
        elif dirfield == self.FIELD_Z:
            vdir=np.array([0.,0.,1.])
            for z in range (self.N):
                vmag+=incmag
                if bInvert:
                    vv = vdir*(Bmaxmag-vmag)
                else:
                    vv = vdir*vmag
                for x in range(self.N):                    
                    for y in range(self.N):                        
                        self.B[x][y][z]= vv
        
        self._makeMainArrow()
        


    
##    def getArrow(self):
##        if self.dirfield == self.FIELD_X:
##            return np.array(np.array([1.,0.,0.]))
##        elif self.dirfield == self.FIELD_Y:
##            return np.array(np.array([0.,1.,0.]))
##        elif self.dirfield == self.FIELD_Z:
##            return np.array(np.array([0.,0.,1.]))
##        
####        N = self.N
#####        Bvec  =  np.copy(self.B[N//2][N//2][N//2])                
####        BvecUp =self.B[N//2][N//2][0]
####        BvecDown =self.B[N//2][N//2][N-1]
####        Bvec=BvecDown - BvecUp
####        print("Bvec down-up in z",Bvec)
####        Bvec  /=  np.linalg.norm(Bvec)
##        return Bvec
            


            
