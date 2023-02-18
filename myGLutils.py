from OpenGL.GL import *
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from OpenGL.arrays import vbo
import numpy as np

from math import atan2,sqrt
from OpenGL.GLU import *       # OpenGL Utility Library, extends OpenGL functionality

def drawArrow(x1,y1,z1,x2,y2,z2,D,c):

    RADPERDEG = np.pi/180.0
    x=x2-x1
    y=y2-y1
    z=z2-z1
    L=sqrt(x*x+y*y+z*z)
    
    glPushMatrix ()
    glTranslate(x1,y1,z1)
    if x!=0. or y!=0.:        
        glRotate(atan2(y,x)/RADPERDEG,0.,0.,1.)
        glRotate(atan2(sqrt(x*x+y*y),z)/RADPERDEG,0.,1.,0.)
    elif z<0:
        glRotate(180,1.,0.,0.)      
    glColor(c[0], c[1], c[2])
    glTranslate(0,0,L-4*D)

    quadObj = gluNewQuadric ()
    gluQuadricDrawStyle (quadObj, GLU_FILL)
    gluQuadricNormals (quadObj, GLU_SMOOTH)
    gluCylinder(quadObj, 2*D, 0.0, 4*D, 32, 1)
    gluDeleteQuadric(quadObj)
    quadObj = gluNewQuadric ()
    gluQuadricDrawStyle (quadObj, GLU_FILL)
    gluQuadricNormals (quadObj, GLU_SMOOTH)
    gluDisk(quadObj, 0.0, 2*D, 32, 1)
    gluDeleteQuadric(quadObj)
    glTranslate(0,0,-L+4*D)
    quadObj = gluNewQuadric ()
    gluQuadricDrawStyle (quadObj, GLU_FILL)
    gluQuadricNormals (quadObj, GLU_SMOOTH)
    gluCylinder(quadObj, D, D, L-4*D, 32, 1)
    gluDeleteQuadric(quadObj)
    quadObj = gluNewQuadric ()
    gluQuadricDrawStyle (quadObj, GLU_FILL)
    gluQuadricNormals (quadObj, GLU_SMOOTH)
    gluDisk(quadObj, 0.0, D, 32, 1)
    gluDeleteQuadric(quadObj)
    
    glPopMatrix ()

def drawVector(pos,v,vlength,headaspect=0.4,color=[1.,0.,0.]):
    
    posend=np.array(pos)+np.array(v)*vlength
    drawArrow(pos[0], pos[1],pos[2],
              posend[0], posend[1],posend[2],headaspect,color)    

def drawCube(x, y, z,sizex,sizey,sizez):
    glPushMatrix()
    glTranslate(-x, -y, -z)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)	
    glBegin(GL_QUADS);
    glColor(1.0, 1.0, 0.0)
    
    # FRONT
    glVertex(-sizex, -sizey, sizez)
    glVertex(sizex, -sizey, sizez)
    glVertex(sizex, sizey, sizez)
    glVertex(-sizex, sizey, sizez)

    # BACK
    glVertex(-sizex, -sizey, -sizez)
    glVertex(-sizex, sizey, -sizez)
    glVertex(sizex, sizey, -sizez)
    glVertex(sizex, -sizey, -sizez)

    

    # LEFT
    glVertex(-sizex, -sizey, sizez)
    glVertex(-sizex, sizey, sizez)
    glVertex(-sizex, sizey, -sizez)
    glVertex(-sizex, -sizey, -sizez)

    # RIGHT
    glVertex(sizex, -sizey, -sizez)
    glVertex(sizex, sizey, -sizez)
    glVertex(sizex, sizey, sizez)
    glVertex(sizex, -sizey, sizez)

    

    # TOP
    glVertex(-sizex, sizey, sizez)
    glVertex(sizex, sizey, sizez)
    glVertex(sizex, sizey, -sizez)
    glVertex(-sizex, sizey, -sizez)

    # BOTTOM
    glVertex(-sizex, -sizey, sizez)
    glVertex(-sizex, -sizey, -sizez)
    glVertex(sizex, -sizey, -sizez)
    glVertex(sizex, -sizey, sizez)

    glEnd()    
    #glTranslate(x, y, z)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
    glPopMatrix()

def arrayToRGBcmap(colors):
    npoints= len(colors)
    colormap = cm.jet
    normalize = mcolors.Normalize(vmin=np.min(colors), vmax=np.max(colors))
    s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
    rgb = s_map.to_rgba(colors).astype(np.float32)    
    return np.delete(rgb, 3,axis=1)  # reome alpha


def drawVBOs(pos_vbo,col_vbo,numpoints):
        glPointSize(4)
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        
        # Position
        pos_vbo.bind()            
        glVertexPointer(3, GL_FLOAT, 0, pos_vbo)
        pos_vbo.unbind()

        # Color    
        col_vbo.bind()
        glColorPointer(3, GL_FLOAT, 0, col_vbo)
        col_vbo.unbind()
                    
        # Draw
        glDrawArrays(GL_POINTS, 0, numpoints)
        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_COLOR_ARRAY)    
def getOglVBOfromArray(array3d,oglscale,minumbral,minpointsogl,maxpointsogl):            
     # select points from array3d over umbral
    
    maxP= np.max(array3d)
    lenArray=np.shape(array3d)
    
    # loop until we get enougth points changing minumbral if required 
    while True:
        umbral=minumbral * maxP
        points =np.where( array3d > umbral)
        numpoints = len(points[0])
        if numpoints<minpointsogl:
            minumbral*=0.9         
        elif numpoints>maxpointsogl:
            minumbral*=1.1
        else:
            break        
    colors=arrayToRGBcmap(array3d[points])
    
    # convert points from tuple to np array
    points = np.swapaxes(points,0,1).astype(np.float32) * oglscale
    
    # swap columns Y Z for ogl cam  
    
    #points[:, [2, 1]] = points[:, [1, 2]]
    
    
    npoints=len(points)
        
    pos_vbo = vbo.VBO(data=points, usage=GL_DYNAMIC_DRAW, target=GL_ARRAY_BUFFER)
    col_vbo = vbo.VBO(data=colors, usage=GL_DYNAMIC_DRAW, target=GL_ARRAY_BUFFER)

    
    
    return pos_vbo,col_vbo,npoints
    
