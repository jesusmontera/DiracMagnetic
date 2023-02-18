

from PyQt5 import QtCore, QtGui, QtWidgets
from dlgdirac import Ui_DlgDiract
import numpy as np

class mydlgDirac(QtWidgets.QDialog):
    """Employee dialog."""
    def __init__(self, parent=None):
        super().__init__(parent)
        # Create an instance of the GUI
        self.ui = Ui_DlgDiract()        
        self.ui.setupUi(self)  # Run the .setupUi() method to show the GUI
        
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
    def getBmodeSpin(self):
        return self.ui.listBspin.currentRow()        
    def getInitialSpin(self):
        if self.ui.opcspin0.isChecked():
            return None
        else:
            cup=complex(float(self.ui.txspinupreal.text()) + float(self.ui.txspinupimag.text()) *1j)
            cdown=complex(float(self.ui.txspindownreal.text()) +float( self.ui.txspindownimag.text()) *1j)
            spin0=np.array( [cup,cdown],np.complex128)
            spin0 = spin0 / np.linalg.norm(spin0)
            
        return spin0
