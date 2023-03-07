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
        self.Bmaxmag=1000.
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
                vpos=np.array([0,0, self.N*0.72])
                vdir=np.array([-1.,0.,0.])
            elif self.dirfield == self.FIELD_Y:
                vpos=np.array([ self.N*0.8,0,0])
                vdir=np.array([0,-1,0])
                
            elif self.dirfield == self.FIELD_Z:
                vpos=np.array([ -self.N*0.8,0,0])
                vdir=np.array([0,0,-1])
        else:
            if self.dirfield == self.FIELD_X:
                vpos=np.array([0,0, self.N*0.72])
                vdir=np.array([1,0,0])
            elif self.dirfield == self.FIELD_Y:
                vpos=np.array([ self.N*0.8,0,0])
                vdir=np.array([0,1,0])
            elif self.dirfield == self.FIELD_Z:
                vpos=np.array([ -self.N*0.8,0,0])
                vdir=np.array([0,0,1.])
        self.mainArrow.append(vpos)
        self.mainArrow.append(vdir)
                
            
            
    def makeB(self,dirfield,Bmaxmag,bInvert):
        self.btext=""
        self.B = np.zeros((self.N,self.N,self.N,3))        
        self.dirfield=dirfield
        self.Bmaxmag=Bmaxmag
        self.bInvert=bInvert
        if not bInvert:
            vmag=Bmaxmag            
            incmag=-2.0 * Bmaxmag/self.N
        else:            
            vmag=-Bmaxmag            
            incmag=2.0 * Bmaxmag/self.N
            
        if dirfield == self.FIELD_X:
            if bInvert:
                vdir=np.array([1.,0.,0.])
            else:
                vdir=np.array([-1.,0.,0.])
            for x in range (self.N):
                vmag+=incmag
                vv = vdir*vmag
                for y in range(self.N):                    
                    for z in range(self.N):                        
                        self.B[x][y][z]= vv
        
        elif dirfield == self.FIELD_Y:
            
            if bInvert:
                vdir=np.array([0.,1.,0.])
            else:
                vdir=np.array([0.,-1.,0.])
                
            for y in range (self.N):
                vmag+=incmag
                vv = vdir*vmag
                for x in range(self.N):                    
                    for z in range(self.N):                        
                        self.B[x][y][z]= vv
        elif dirfield == self.FIELD_Z:
            if bInvert:
                vdir=np.array([0.,0.,1])
            else:
                vdir=np.array([0.,0.,-1])
            for z in range (self.N):
                vmag+=incmag                
                vv = vdir*vmag
                for x in range(self.N):                    
                    for y in range(self.N):                        
                        self.B[x][y][z]= vv
        
        self._makeMainArrow()
        


    

            
