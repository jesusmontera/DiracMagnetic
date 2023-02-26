

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
from dlgdirac import Ui_DlgDiract
from GLBlochWidget import GLBlochWidget
from auxfunctions import getBlochVector



class mydlgDirac(QtWidgets.QDialog):
    """Employee dialog."""
    def __init__(self, parent=None):
        super().__init__(parent)
        # Create an instance of the GUI
        self.ui = Ui_DlgDiract()        
        self.ui.setupUi(self)  # Run the .setupUi() method to show the GUI
        self.ui.btPlotSpin.clicked.connect(self.PlotOGLSpin)
        self.glWidget = GLBlochWidget()
        self.glWidget.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.ui.oglblochlayout.addWidget(self.glWidget)        
        #self.glWidget.setspacesize( 4)
        self.glWidget.camgoto([0,0, - 25.])
        self.PlotOGLSpin()
        
    def PlotOGLSpin(self):
        
        spin=self.getInitialSpin()
        sdir=getBlochVector(spin)
        self.glWidget.clearArrows()
        spos=np.array([5.,5.,5.])
        scolor=np.array([0.,0.8,0.])
        self.glWidget.addArrow(spos,sdir,scolor)
        self.glWidget.update()
        
        
            
            
    def myvalidator(self,w,vmin,vmax,vdefault=0.0):        
        s=w.text()
        if s=="":
            w.setText(str(vdefault))
        else:
            v=  float(s)
            if v < vmin or v > vmax:
                QtWidgets.QMessageBox.about(self, "Value incorrect",
                                            "limit for " +w.objectName() + " is from " + str(vmin) + " to " + str(vmax))
                return False        
        return True
    
        
    def accept(self):
        if not self.myvalidator(self.ui.txposX,-0.8,0.8,0.0): return
        if not self.myvalidator(self.ui.txposY,-0.8,0.8,0.0): return
        if not self.myvalidator(self.ui.txposZ,-0.8,0.8,0.0): return

        if not self.myvalidator(self.ui.txkX,-137,137): return
        if not self.myvalidator(self.ui.txkY,-137,137): return
        if not self.myvalidator(self.ui.txkZ,-137,137): return
        
        super().accept()     
        
        
    def getPos(self):
        x = float( self.ui.txposX.text())
        y = float( self.ui.txposY.text())
        z = float( self.ui.txposZ.text())
        return np.array([x,y,z])
    def getK(self):
        x = float( self.ui.txkX.text())
        y = float( self.ui.txkY.text())
        z = float( self.ui.txkZ.text())
        return np.array([x,y,z])
    def getMaxFrames(self):
        return int(self.ui.txMaxSteps.text())
    
    def getInitialSpin(self):
        if self.ui.opcspin0.isChecked():
            return None
        else:
            cup=complex(float(self.ui.txspinupreal.text()) + float(self.ui.txspinupimag.text()) *1j)
            cdown=complex(float(self.ui.txspindownreal.text()) +float( self.ui.txspindownimag.text()) *1j)
            spin0=np.array( [cup,cdown],np.complex128)
            spin0 = spin0 / np.linalg.norm(spin0)
            
        return spin0
