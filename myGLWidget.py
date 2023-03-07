
from PyQt5 import QtGui
from PyQt5.QtCore import Qt
from PyQt5 import QtOpenGL    # provides QGLWidget, a special OpenGL QWidget
from OpenGL.GL import *        # python wrapping of OpenGL
from OpenGL import GLU        # OpenGL Utility Library, extends OpenGL functionality
from OpenGL.arrays import vbo
import numpy as np
from myGLutils import drawCube,drawVector,drawVBOs

class GLWidget(QtOpenGL.QGLWidget):    
    def __init__(self, parent=None):        
        self.parent = parent        
        self.pos_vbo=None
        self.col_vbo=None
        self.numpoints=0
        self.spacesize = 64.0
        self.campos=np.array([0, 0, 0],np.float32)
        self.camYangle=0.0
        self.labelinfo=None
        self.arrows=[]
        self.Btext=""
        self.bTextInfo=True
        # camera speed for rotating and moving
        self.rotangle=3.0
        self.speed=2.0
        QtOpenGL.QGLWidget.__init__(self, self.parent)
    def setTextInfo(self,bShow):
        self.bTextInfo=bShow
        
    def drawOGLText(self,pos,text,fontsize=12, color=[1.,1.,1.]):
        glPushMatrix ()        
        glTranslate(-pos[0], -pos[1], -pos[2])
        glColor3f(color[0], color[1],color[2])
        font = QtGui.QFont("Arial");
        font.setPointSize(fontsize)            
        self.renderText(0,0,0,text,font)
        glPopMatrix ()
        
    def scaleArrow(self,na,s):
        self.arrows[na*4+3]*=s
    
    def clearArrows(self):
        self.arrows=[]    
    def setspacesize(self,spacesize):
        self.spacesize = spacesize
        self.update()
        print("self.spacesize",self.spacesize)
    def setVBOs(self,pos_vbo,col_vbo,numpoints):
        self.pos_vbo=pos_vbo
        self.col_vbo=col_vbo
        self.numpoints=numpoints
        #self.update()            
        
    def addArrow(self,apos, avector, color,length=22.):        
        self.arrows.append(apos)
        self.arrows.append(avector)
        self.arrows.append(color)
        self.arrows.append(length)

        
        
    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Up:
            if event.modifiers() == Qt.ShiftModifier:
                self.movecam(Yinc=-self.speed)
            else:
                self.movecam(0,self.speed)                
                                                
        elif event.key() == Qt.Key_Down:
            if event.modifiers() == Qt.ShiftModifier:
                self.movecam(Yinc=self.speed)                
            else:
                self.movecam(0,-self.speed)
                
                
        elif event.key() == Qt.Key_Right:
            if event.modifiers() == Qt.ShiftModifier:
                self.movecam(Ry=self.rotangle)
            else:
                self.movecam(self.speed)                
                
                                
        elif event.key() == Qt.Key_Left:
            if event.modifiers() == Qt.ShiftModifier:
                self.movecam(Ry=-self.rotangle)
            else:
                self.movecam(-self.speed)                
                
                
                
                
    def setlabelinfo(self, labelinfo):
        self.labelinfo=labelinfo
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
        if self.labelinfo is not None:            
            self.labelinfo.setText(str(np.round(self.campos,1)))

    def drawArrows(self):
        for i in range(len( self.arrows)//4):
            self.drawArrow(i)
            
    def drawArrow(self,i):
        
        if len( self.arrows)//4 > i:                                    
            pos=self.arrows[i*4]
            vec=self.arrows[i*4+1]
            color=self.arrows[i*4+2]
            length=self.arrows[i*4+3]
            head= length*0.025
            drawVector(pos,vec,length,head,color)                    
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
        # Y
        
        
    def paintGL(self):
        
        #self.makeCurrent()
                
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)                
        glLoadIdentity()
        # convert ogl axis(Y up x right) to matplotlib axis (Z up y right)
        # rotating 
        #glRotate(-180., 0.0, 0.0, 1.0)
        glRotate(90+self.camYangle, 0.0, 1.0, 0.0)            
        
        glTranslate(self.campos[0], self.campos[1], self.campos[2])
        glRotate(90, 1.0, 0.0, 0.0)            
                
        drawCube(0,0, 0,self.spacesize/2,self.spacesize/2,self.spacesize/2)

        self.drawAxes()
        if self.Btext!="" and self.bTextInfo:            
            pos = [ 0,-self.spacesize*0.7,0]
            self.drawOGLText(pos,self.Btext,fontsize=14,color=[0.9, 0, 0])
       
        if self.pos_vbo is not None:

            # B arrow
            if len( self.arrows)//4==2:
                self.drawArrow(0)
            # translate to psi coordinates 
            glTranslate(-self.spacesize/2.,
                        -self.spacesize/2., -self.spacesize/2)             
            # spin arrow 
            if len( self.arrows)//4==2:
                self.drawArrow(1) 
            else:  self.drawArrow(0)

            drawVBOs( self.pos_vbo, self.col_vbo , self.numpoints)
            
            
            ########################        
        else:
            if self.bTextInfo:
                self.drawOGLText([0,self.spacesize/2,0],"OPENGL DIRAC MAGNETIC FIELD")
                self.drawOGLText([0,self.spacesize/2,10],"Click and move with arrow keys + shift")            
                                
            self.drawArrows()
            

        
        
        
        #glPopMatrix()
        
        
        #painter.endNativePainting()
##        painter = QtGui.QPainter(self)
##        painter.begin(self)
##        pen = QtGui.QPen(QtGui.QColor(255,0,0), 10)
##        painter.setPen(pen)        
##        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
##        painter.setRenderHint(QtGui.QPainter.TextAntialiasing,True)
##        painter.drawText(10,10, "texto QPainter en QGLWidget.PaintGL")
##        painter.end()
  
    def camgoto(self,campos):
        self.campos=np.array(campos,np.float32)
        self.update()
        if self.labelinfo is not None:            
            self.labelinfo.setText(str(np.round(self.campos,1)))
    
