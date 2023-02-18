

from PyQt5 import QtWidgets
from dlgogl import Ui_DlgOgl

class mydlgOgl(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super().__init__(parent)
        # Create an instance of the GUI
        self.ui = Ui_DlgOgl()
        # Run the .setupUi() method to show the GUI        
        self.ui.setupUi(self)            
        
    
