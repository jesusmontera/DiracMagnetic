

from PyQt5 import QtCore, QtGui, QtWidgets
from dlgschrod import Ui_DlgSch
import numpy as np

class mydlgSch(QtWidgets.QDialog):
    """Employee dialog."""
    def __init__(self, parent=None):
        super().__init__(parent)
        # Create an instance of the GUI
        self.ui = Ui_DlgSch()
        # Run the .setupUi() method to show the GUI
        
        self.ui.setupUi(self)
        
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
        

        if not self.myvalidator ( self.ui.txSpinX, -1. , 1. ): return
        if not self.myvalidator ( self.ui.txSpinY, -1. , 1. ): return
        if not self.myvalidator ( self.ui.txSpinZ, -1. , 1. ): return
        
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
        blochspin= np.zeros(3)
        blochspin[0]=float(self.ui.txSpinX.text())
        blochspin[1]=float(self.ui.txSpinY.text())
        blochspin[2]=float(self.ui.txSpinZ.text())
        norm = np.linalg.norm(blochspin)        
        return blochspin/norm    
        
        
