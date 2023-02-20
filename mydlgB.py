# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Bdlg.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from dlgB import Ui_Dialog
from myMagnetic import magneticlass
import numpy as np
from GLBWidget import GLBWidget
class Bdlg(QtWidgets.QDialog):
    """Employee dialog."""
    def __init__(self, parent):
        super().__init__(parent)
        # Create an instance of the GUI
        self.ui = Ui_Dialog()
        # Run the .setupUi() method to show the GUI        
        self.ui.setupUi(self)
        #########  menu #######################3
        self.parent=parent
        bar = QtWidgets.QMenuBar()
        file = bar.addMenu("File")
        file.addAction("Open B file",self.openBfile)
        file.addAction("Save B file",self.saveBfile)
        self.ui.menuLayout.setMenuBar(bar)
        
        
        #########  B objet and matplolib  #######################3
        
        self.ui.btMakeB.clicked.connect(self.makeB)
                        
        self.glWidget = GLBWidget()
        self.glWidget.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.ui.oglayout.addWidget(self.glWidget)        
        self.glWidget.setspacesize( parent.N * parent.oglscale)
        self.glWidget.camgoto([0,0, - parent.N * parent.oglscale * 2])
        
    def openBfile(self):        
        file , check = QtWidgets.QFileDialog.getOpenFileName(self,
                                                             "Open B npy file",
                                                             "", "npy Files (*.npy)")
        if check:            
            #self.parent.B.B = np.load(file,allow_pickle=True)
            self.parent.B.load_npyFile(file)
            self.plotB()
    def saveBfile(self):
        if self.parent.B.B is None: return
        file , check = QtWidgets.QFileDialog.getSaveFileName(self,
                                                  'save B npy file',
                                                  "", "npy Files (*.npy)")
        if check:
            np.save(file, self.parent.B.B)            
    def getBiggestMag():        
        np.linalg.norm(v)
    def plotB(self):                
        factoreduce=10 # reduce number of arrows so we can see then ok
        self.glWidget.clearArrows()
        color=np.array([0.8,0.,0.])
        N=self.parent.N

        # this 3 lines are for arrows  scale after adding them
        maxarrows=int((N/factoreduce+1)**3)
        magnitudes=np.zeros(maxarrows)
        numarrows=0
        # add the arrows to glwidget
        for x in range (N):
            if x % factoreduce ==0:
                for y in range(N):
                    if y % factoreduce ==0:
                        for z in range(N):
                            if z % factoreduce ==0:
                                vdir=np.copy(self.parent.B.B[x][y][z])                                
                                norm=np.linalg.norm(vdir)
                                
                                if self.parent.B.bInvert:
                                    vdir*=-1
                                    
                                if norm > 0:                                                                        
                                    pos=np.array([x,y,z])
                                    self.glWidget.addArrow(pos, vdir/norm, color)
                                    magnitudes[numarrows]=norm
                                    numarrows+=1
        #now scale the arrows
        maxmag= np.max(magnitudes)
        for i in range(numarrows):
            s=magnitudes[i]/maxmag
            self.glWidget.scaleArrow(i,s*9)

        minmag= np.round(np.min(magnitudes),3)
        maxmag=np.round(maxmag,3)
        self.ui.lbinfo.setText("B magnitude from " + str(minmag) +  " to " + str(maxmag) + " " + self.parent.B.btext)
        
        self.glWidget.update()
                                    
    
    def makeB(self):                
        Baxis= self.ui.listdir.currentRow()
        Bmaxmag=float(self.ui.txmaxB.text())
        self.parent.B.clear()
        self.ui.lbinfo.setText("Making magnetic B array of vectors")
        
        bInvert=self.ui.checkInvert.isChecked()
        self.parent.B.makeB(Baxis,Bmaxmag,bInvert)                
        self.plotB()
        
        
    
