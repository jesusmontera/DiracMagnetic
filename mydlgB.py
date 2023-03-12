
from PyQt5 import QtCore, QtGui, QtWidgets
from dlgB import Ui_Dialog
from myMagnetic import magneticlass
import numpy as np
from myGLWidget import GLWidget
import matplotlib.pyplot as plt
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
                        
        self.glWidget = GLWidget()
        self.glWidget.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.ui.oglayout.addWidget(self.glWidget)        
        self.glWidget.setspacesize( parent.N * parent.oglscale)
        self.glWidget.camgoto([0,0,-parent.N * parent.oglscale * 2.])
        
        
        
    def openBfile(self):        
        file , check = QtWidgets.QFileDialog.getOpenFileName(self,
                                                             "Open B npy file",
                                                             "", "npy Files (*.npy)")
        if check:                        
            self.parent.B.load_npyFile(file)            
            #self.parent.B.B *=4.
##            self.RotateB()
            self.plotB(True)
    def saveBfile(self):
        if self.parent.B.B is None: return
        file , check = QtWidgets.QFileDialog.getSaveFileName(self,
                                                  'save B npy file',
                                                  "", "npy Files (*.npy)")
        if check:
            np.save(file, self.parent.B.B)
    
        
    def getBiggestMag():        
        np.linalg.norm(v)
    def RotateB(self):
        if self.parent.B.B is not None:            
            self.parent.B.B = np.rot90(self.parent.B.B, k = 2, axes = (0, 1))
            self.plotB()
    def plotB(self, bfromnpy=False):                
        factoreduce=10 # reduce number of arrows so we can see then ok
        self.glWidget.clearArrows()
        color=np.array([0.8,0.,0.])
        N=self.parent.N

        # this 3 lines are for arrows  scale after adding them
        maxarrows=int((N/factoreduce+1)**3)
        magnitudes=np.zeros(maxarrows)
        numarrows=0
        # add the arrows to glwidget
        B = self.parent.B.B
        wc=self.glWidget.spacesize/2.
        for x in range (N):
            if x % factoreduce ==0:
                for y in range(N):
                    if y % factoreduce ==0:
                        for z in range(N):
                            if z % factoreduce ==0:                                
                                vdir=np.copy(B[x][y][z])  
                                norm=np.linalg.norm(vdir)                                                                
                                if norm > 0:
                                    #vdir*=-1.    
                                    pos=np.array([x-wc,y-wc,z-wc])
                                    self.glWidget.addArrow(pos, vdir/norm, color,1.)
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
        if bfromnpy:
            self.ui.txmaxB.setText(str(maxmag))
        self.glWidget.update()
                                    
    
    def makeB(self):
        sResult=""
        Bmaxmag=float(self.ui.txmaxB.text())
        self.ui.btMakeB.setText("wait")
        self.ui.lbinfo.setText("Making magnetic B array of vectors...wait")        
        QtGui.QGuiApplication.processEvents()
        if self.ui.opcBApp.isChecked():            
            Baxis= self.ui.listdir.currentRow()            
            self.parent.B.clear()            
            self.parent.B.makeB(Baxis,Bmaxmag)            
            self.plotB()
            sResult="B maked ok "
        else:
            try:
                plt.close()
                from BmagpyGerlachZ import MakeBmagpyZ
                L_au=float(self.parent.txL.text())
                self.parent.B.B = MakeBmagpyZ(N = self.parent.N,                                                   
                                                   L_au= L_au,
                                                   B0 = Bmaxmag, Bdir=-1)
                
                sResult="B magpy non uniform field in Z maked ok "
                self.plotB()
                self.parent.B.mainArrow=[]
                self.parent.B.btext="magpy"
                
            except ImportError:
                print("ERROR :exception becuase u don't have magpylib")
                print("solution: load the npy sample file BgerlachZ.npy with the menu")
                sResult="error: try loading BgerlachZ.npy with the menu"
                pass
        self.ui.btMakeB.setText("Make B")
        self.ui.lbinfo.setText(sResult)
        print("end making B")
            
        
        
        
    
