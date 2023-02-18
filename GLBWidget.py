
from PyQt5 import QtGui
from PyQt5.QtCore import Qt
from PyQt5 import QtOpenGL    # provides QGLWidget, a special OpenGL QWidget
from OpenGL.GL import *        # python wrapping of OpenGL
from OpenGL import GLU        # OpenGL Utility Library, extends OpenGL functionality
from OpenGL.arrays import vbo
import numpy as np
from myGLutils import drawCube,drawVector

class GLBWidget(QtOpenGL.QGLWidget):    
    def __init__(self, parent=None):        
        self.parent = parent                
        self.numpoints=0
        self.spacesize = 64.0
        self.campos=np.array([0, 0, 0],np.float32)
        self.camYangle=0.0
        
        self.arrows=[]
        
        # camera speed for rotating and moving
        self.rotangle=3.0
        self.speed=2.0
        QtOpenGL.QGLWidget.__init__(self, self.parent)
                                    
    def setspacesize(self,spacesize):
        self.spacesize = spacesize
        self.update()
        print("self.spacesize",self.spacesize)
    def clearArrows(self):
        self.arrows=[]
    def addArrow(self,apos, avector, color):        
        self.arrows.append(apos-self.spacesize/2)
        self.arrows.append(avector)
        self.arrows.append(color)

        
        
    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Up:
            if event.modifiers() == Qt.ShiftModifier:                
                self.movecam(Yinc=-self.speed)
            else:
                self.movecam(-self.speed)                                
        elif event.key() == Qt.Key_Down:
            if event.modifiers() == Qt.ShiftModifier:                
                self.movecam(Yinc=self.speed)
            else:
                self.movecam(self.speed)
                
        elif event.key() == Qt.Key_Right:
            if event.modifiers() == Qt.ShiftModifier:
                self.movecam(Ry=self.rotangle)
            else:
                self.movecam(0,self.speed)                
        elif event.key() == Qt.Key_Left:
            if event.modifiers() == Qt.ShiftModifier:
                self.movecam(Ry=-self.rotangle)
            else:
                self.movecam(0,-self.speed)
                    
    def __del__(self):
        print("GLWidget destructor")
        
        
    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(0., 0, 0))    # initialize the screen to blue        
        glEnable(GL_DEPTH_TEST)                  # enable depth testing        
        
        
         
    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        aspect = width / float(height)

        GLU.gluPerspective(45.0, aspect, 1.0, 300.0)
        glMatrixMode(GL_MODELVIEW)
    def movecam(self,speedfrontal=0,speedlateral=0,Yinc=0, Ry=0):
        
        if Yinc != 0.: self.campos[1] +=Yinc        
        if Ry  != 0.:
            self.camYangle+=Ry
            if self.camYangle>=360: self.camYangle-=360.
            if self.camYangle<0: self.camYangle+=360.
            
        # forward move
        dirx= np.sin( self.camYangle * np.pi/180.0)
        dirz= -np.cos(self.camYangle*np.pi/180.0)
        self.campos[0]  += dirx * speedfrontal
        self.campos[2]  += dirz * speedfrontal
        # lateral move
        angulolateral =self.camYangle-90.0
        if angulolateral<0:
            angulolateral+=360
            
        dirx= np.sin(angulolateral*np.pi/180.0)
        dirz= -np.cos(angulolateral*np.pi/180.0)
        self.campos[0]  +=dirx * speedlateral
        self.campos[2]  +=dirz * speedlateral	
        
        #self.paintGL()
        self.update()
        
            
    def drawArrows(self):        
        na=len(self.arrows)//3        
        for i in range (na):        
            pos=self.arrows[i*3]
            vec=self.arrows[i*3+1]
            color=self.arrows[i*3+2]                            
            drawVector(pos,vec,5.,0.3,color)                    
    def drawAxes(self):        
        color=[1.,1.,1.]
        # X
        pos = [self.spacesize/2, self.spacesize/2+6, -self.spacesize/2]
        self.drawOGLText(pos,"X ",fontsize=16)
        pos =np.array(pos)*-1.
        pos[0]+=3.5
        pos[1]+=1.3
        vec=[1,0.,0]
        drawVector(pos,vec,7,0.3,color)
        # Y
        pos = [-self.spacesize/2-4,self.spacesize/2,-self.spacesize/2]
        self.drawOGLText(pos,"Y",fontsize=16)
        pos =np.array(pos)*-1.
        pos[1]+=4.
        pos[0]+=1.
        vec=[0,1.0,0]
        drawVector(pos,vec,7,0.3,color)
    def drawOGLText(self,pos,text,fontsize=12):
        glPushMatrix ()        
        glTranslate(-pos[0], -pos[1], -pos[2])
        glColor3f(1., 1.,1.)
        font = QtGui.QFont("Arial");
        font.setPointSize(fontsize)            
        self.renderText(0,0,0,text,font)
        glPopMatrix ()    
    def paintGL(self):
        
        #self.makeCurrent()
                
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)                
        glLoadIdentity()
        glRotate(self.camYangle, 0.0, 1.0, 0.0)            
        glTranslate(self.campos[0], self.campos[1], self.campos[2])
                
        drawCube(0,0, 0,self.spacesize/2,self.spacesize/2,self.spacesize/2)
        self.drawAxes()
        self.drawArrows()
        
  
    def camgoto(self,campos):
        self.campos=np.array(campos,np.float32)
        self.update()
        
    
