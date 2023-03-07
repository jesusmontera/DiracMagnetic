
import numpy as np
from time import time
from threading import Thread
from PyQt5 import QtCore, QtGui, QtWidgets
from mainwnd import Ui_MainWindow
from myGLWidget import GLWidget
from myGLutils import getOglVBOfromArray
from mydlgogl import mydlgOgl

from myDirac3D import myDirac3D # encapsulates dirac (both split step and free eigen periodic)
from myPauli3D import myPauli3D # encapsulates Pauli class
from mydlgdirac import mydlgDirac
from myMagnetic import magneticlass
from mydlgB import Bdlg
from mydlgSch import mydlgSch

from myScroduinger3D import mySchroduinger3D

class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent=parent)
        self.setupUi(self)
        self.N=int(self.txN.text())        
        self.oglscale=1.0
        self.busy=False
        self.playType=""
        
        ############ add OpenGL widget to oglalyout ######################333
        self.glWidget = GLWidget()
        self.glWidget.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.oglalyout.addWidget(self.glWidget)
        self.glWidget.setlabelinfo(self.lbinfo)
        self.glWidget.setspacesize( self.N * self.oglscale)
        
        self.glWidget.camgoto([0,0,-self.N * self.oglscale * 2.])
        
        #### for extracting  N OpenGL points from Dirac probability numpy array
        ### can be changed in OGL dialog launched with Play OGL
        self.minumbralogl=0.0002
        self.minpointsogl=800
        self.maxpointsogl=3000
        
        ############  Dirac and schroduinger objects ######################333
        self.Dirac=myDirac3D(self.N)
        self.Pauli=myPauli3D(self.N)
        self.schroduinger= mySchroduinger3D(self.N)
        self.btMake.clicked.connect(self.makePSI)
        self.btPlay.clicked.connect(self.btPlayOGL)        
        #################################################        
        self.B = magneticlass(self.N) # magnetic object
        #################################################
        self.sliderOGL.setRange(0, 0)                
        self.mnuSaveDirac.triggered.connect(self.saveDiracFile)
        self.mnuLoadDirac.triggered.connect(self.loadDiracFile)
        self.mnuSaveSch.triggered.connect(self.saveSchFile)
        self.mnuLoadSch.triggered.connect(self.loadSchFile)
        
        self.mnumagshow.triggered.connect(self.showDlgB)
        self.sliderOGL.sliderMoved['int'].connect(self.renderOGLframe)
        self.dlgdirac= None
        self.dlgmag= None
        self.dlgscho= None
        
        
        #################################################
                    
            
    def showDlgB(self):
        # magnetic file dialog                       
        if self.dlgmag is None:            
            self.dlgmag= Bdlg(self)
        Bbefore = self.B.B
        if self.dlgmag.exec():            
            self.schroduinger.clear()
            self.Dirac.clear()            
            self.btPlay.setEnabled(False)
            self.glWidget.setVBOs(None,None,0)
            self.sliderOGL.setSliderPosition(0)
            self.sliderOGL.setRange(0, 0)
            self.glWidget.arrows=[]            
                # add B arrow to OpenGL
            self.glWidget.Btext=self.B.btext            
            if self.B.mainArrow  != []:                                                           
                self.glWidget.addArrow(self.B.mainArrow[0] , self.B.mainArrow[1],[0.8,0.,0.],11.)
        else:
            self.B.B=Bbefore
            
        self.glWidget.update()
        
                                        
                
    def saveSchFile(self):
        file , check = QtWidgets.QFileDialog.getSaveFileName(self,
                                                  'save Schroduinger npy file',
                                                  "", "npy Files (*.npy)")
        if check:
            self.schroduinger.saveFileProb(file)
            
    def saveDiracFile(self):        
        file , check = QtWidgets.QFileDialog.getSaveFileName(self,
                                                  'save Dirac npy file',
                                                  "", "npy Files (*.npy)")
        if check:
            self.Dirac.saveFileProb(file)
                                

    def loadSchFile(self):
        file , check = QtWidgets.QFileDialog.getOpenFileName(self,
                                                             "Open Schroduinger npy file",
                                                             "", "npy Files (*.npy)")
        if check:            
            self.schroduinger.clear()
            self.Dirac.clear()
            self.B.clear()
            self.radioSchro.setChecked(True)
            animframes = self.schroduinger.loadFileProb(file)        
            self.restartOGLcontrols(animframes)
            print("file loaded ok with ",animframes, " frames")

    def restartOGLcontrols(self,animframes):            
            self.sliderOGL.setRange(0, animframes)
            self.sliderOGL.setSliderPosition(0)
            self.sliderOGL.update()            
            self.lbinfo.setText("Loaded ok. animation pos 0 from " + str(animframes))
            QtGui.QGuiApplication.processEvents()
            self.glWidget.camYangle=0.0
            
            self.glWidget.camgoto([0,0,-self.N * self.oglscale * 2.])
            
            
            self.glWidget.arrows=[]
            self.renderOGLframe()
            self.btPlay.setEnabled(True)
        
    def loadDiracFile(self):
        file , check = QtWidgets.QFileDialog.getOpenFileName(self,
                                                             "Open Dirac npy file",
                                                             "", "npy Files (*.npy)")
        if check:            
            self.schroduinger.clear()
            self.Dirac.clear()
            self.B.clear()
            self.radioDirac.setChecked(True)
            animframes = self.Dirac.loadFileProb(file)
            self.restartOGLcontrols(animframes)
            print("file loaded ok with ",animframes, " frames")
                                    
        
    def renderOGLframe(self):
                
        animpos = self.sliderOGL.sliderPosition()
        animframes = self.sliderOGL.maximum()        
        
        if animpos < animframes:
            # probability
            if self.radioDirac.isChecked():
                obj = self.Dirac
                ss="dirac "
            elif self.radioPauli.isChecked():
                obj = self.Pauli
                ss="Pauli "
            else:
                obj = self.schroduinger
                ss="schroduinger "

            prob = obj.getProbability(animpos)
            
            points, colors,npoints = getOglVBOfromArray(array3d= prob,
                                                        oglscale= self.oglscale,
                                                        minumbral=self.minumbralogl,
                                                        minpointsogl=self.minpointsogl,
                                                        maxpointsogl=self.maxpointsogl)
            
            self.glWidget.setVBOs(points, colors,npoints)
            # arrows from spin and B magnetic
            # spin from actual dirac wf

            self.glWidget.arrows=[]

            spinpos, spinvec = obj.getSpinbloch(animpos) 
            
            if self.B.mainArrow != []:
                self.glWidget.addArrow(self.B.mainArrow[0], self.B.mainArrow[1],[0.8,0.,0.],11.)                             
                
            self.glWidget.addArrow(spinpos, spinvec,[0.,0.7,0.])
            # paint open GL
            self.glWidget.update()
            # update label widgets in main window            
            self.lbinfo.setText(ss + str(animpos) + " from " + str(animframes) +
                                " spin ="  + str(np.round(spinvec,3)) +
                                "\nOGL points =" + str(len(points)))
            return True
        else:
            return False
    
    def threadPlayOGL(self):
        self.sliderOGL.setSliderPosition(0)                
        self.btPlay.setText("Stop")
        self.btMake.setEnabled(False)
        QtGui.QGuiApplication.processEvents()                
        while self.busy:
            tstart = time()
            if not self.renderOGLframe(): break            
            animpos = self.sliderOGL.sliderPosition()
            self.sliderOGL.setSliderPosition(animpos+1)
            self.sliderOGL.update()
            QtGui.QGuiApplication.processEvents()
            while time()-tstart < 0.1: # 10 fotos per second
                if not self.busy: break
                QtGui.QGuiApplication.processEvents()
            
        self.sliderOGL.setSliderPosition(0)        
        self.btPlay.setText("Play OpenGL")
        self.btMake.setEnabled(True)
        self.busy=False
        QtGui.QGuiApplication.processEvents()                                                
        
    def threadmakePSI(self,obj,maxframes):
        #obj can be Dirac Pauli or Schroduinger object
        self.animframes=0
        self.sliderOGL.setSliderPosition(0)
        self.btPlay.setEnabled(False)        
        for i in range(maxframes):
            if not self.busy: break
            tstart = time()
            obj.doAnimFrame(Bmagnetic=self.B.B )
            tframe = time()-tstart
            self.lbinfo.setText("frame " + str(i+1) + " took " + str(round(tframe,2)) + " seconds")
            QtGui.QGuiApplication.processEvents()                                                        
        self.sliderOGL.setRange(0, i)
        self.btPlay.setEnabled(True)
        self.btMake.setText("Make")        
        self.busy=False
        if i >0:
            self.renderOGLframe()
            self.btPlay.setEnabled(True)
    def makePSI(self):        
        
        self.busy= not self.busy
        if not self.busy: return
        
        i=0
        if self.radioDirac.isChecked() or  self.radioPauli.isChecked():
            # Dirac and Pauli use the same dialog box but Pauli don't have
            # split method so we'll hide the group box widget that choose split or eigen
            
            if self.dlgdirac is None:
                self.dlgdirac = mydlgDirac(self)                            

            if self.radioPauli.isChecked():                                    
                self.dlgdirac.ui.groupBoxMode.setVisible(False)
                self.dlgdirac.ui.groupBoxB.setVisible(True)
                if self.B.B is None:
                    self.dlgdirac.setWindowTitle("Pauli options (B magnetic isn't loaded)")                
                    self.dlgdirac.ui.groupBoxB.setEnabled(False)
                else:
                    self.dlgdirac.setWindowTitle("Pauli options (B magnetic is loaded)")                
                    self.dlgdirac.ui.groupBoxB.setEnabled(True)
                    
            else:
                self.dlgdirac.setWindowTitle("Dirac options")
                self.dlgdirac.ui.groupBoxMode.setVisible(True)
                self.dlgdirac.ui.groupBoxB.setVisible(False)
            if self.dlgdirac.exec():
                self.schroduinger.clear()
                self.Pauli.clear()                        
                self.Dirac.clear()                        
                POS0=self.dlgdirac.getPos()
                K0 = self.dlgdirac.getK()
                DT=float(self.txDT.text())
                L = float(self.txL.text())                
                
                spin0 = self.dlgdirac.getInitialSpin()
                
                maxframes=self.dlgdirac.getMaxFrames()
                self.btMake.setText("wait...")
                self.lbinfo.setText("initializing Dirac, wait a few seconds")
                QtGui.QGuiApplication.processEvents()

                if self.radioDirac.isChecked():
                    bSplit=self.dlgdirac.ui.opcSplit.isChecked()                
                    self.Dirac.initAnimation(bSplit,L,DT,K0, POS0, initial_spin=spin0, B=self.B.B)
                    thmakepsi = Thread(target=self.threadmakePSI, args=(self.Dirac,maxframes,))
                else:
                    if self.dlgdirac.ui.opcBpot.isChecked():
                        Bmode="pot"
                    else:
                        Bmode="dot"
                    self.Pauli.initAnimation(L,DT,K0, POS0, initial_spin=spin0,
                                             B = self.B.B , Bmode= Bmode)
                    thmakepsi = Thread(target=self.threadmakePSI, args=(self.Pauli,maxframes,))
                    
                self.btMake.setText("Stop")
                QtGui.QGuiApplication.processEvents()
                
                ### thread                                
                thmakepsi.start()
            else: # calcel dialog button
                self.busy=False
            
        else: # schroduinger dialog box
            if self.dlgscho is None:
                self.dlgscho = mydlgSch(self)                                        
            
            if self.dlgscho.exec():                
                self.schroduinger.clear()
                self.Pauli.clear()                        
                self.Dirac.clear()      
                L = float(self.txL.text())
                DT=float(self.txDT.text())
                POS0=self.dlgscho.getPos()
                K0 = self.dlgscho.getK()                
                
                spinbloch = self.dlgscho.getInitialSpin()
                maxframes=self.dlgscho.getMaxFrames()            
                self.lbinfo.setText("initializing schroduinger")
                self.btMake.setText("Stop")
                QtGui.QGuiApplication.processEvents()
                # shcroduinher in meters L=1e-8,DT=5e-16                
                self.schroduinger.initAnimation(L,DT,K0, POS0, initial_spin=spinbloch,B=self.B.B)
                thmakepsi = Thread(target=self.threadmakePSI, args=(self.schroduinger,maxframes,))
                thmakepsi.start()
            else: # calcel dialog button
                self.busy=False
                            
        
    def btPlayOGL(self):
        self.busy= not self.busy
        if not self.busy: return
        
        dlgogl= mydlgOgl(self)            
        if dlgogl.exec():
            self.minumbralogl=float(dlgogl.ui.txumbral.text())
            self.minpointsogl=int(dlgogl.ui.txminPoints.text())
            self.maxpointsogl=int(dlgogl.ui.txmaxpoints.text())                            
            thplayopgl = Thread(target=self.threadPlayOGL)
            thplayopgl.start()
        
        

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())

