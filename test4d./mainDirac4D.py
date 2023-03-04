import numpy as np
from time import time,sleep
from threading import Thread
from PyQt5 import QtCore, QtGui, QtWidgets
from mainwnd import Ui_MainWindow
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from dirac4d import Dirac4D
from auxfunctions import getBlochVector
from potentials import makeYBarrier
class PlotCanvas(FigureCanvas):
    def __init__(self, parent = None, width = 5, height = 5, dpi = 100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(1,1,1)
 
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        self.figure.tight_layout()
        self.im=None
 
 
    def plot(self,p):            
        if self.im is None:
            self.im = self.axes.imshow(p,
                                   interpolation='bilinear',
                                   origin='lower', cmap='hot')
        else:                    
            self.im.set_data(p)
        self.draw()
 
class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent=parent)
        self.setupUi(self)
        
        self.busy=False
        self.btPlay.clicked.connect(self.funcPlay)
        self.canvasplot = PlotCanvas(self.lbogl, width=4.5, height=4.5)
        self.canvasplot.move(0,0)
        self.dirac=None        
        self.V=None
    def funcPlay(self):
        self.busy=not self.busy        
        QtGui.QGuiApplication.processEvents()
        
        if not self.busy: return
        
        N=int(self.txN.text())
        L=float(self.txL.text())
        self.V = makeYBarrier(N,L*5.29177210903E-11)        
        self.dirac= Dirac4D(N,self.V)        
                
        pos1=np.array([float(self.txX1.text()),float(self.txY1.text())])
        pos2=np.array([float(self.txX2.text()),float(self.txY2.text())])
                       
        k1=np.array([float(self.txKX1.text()),float(self.txKY1.text()),0.])
        k2 = np.array([float(self.txKX2.text()),float(self.txKY2.text()),0.])
        spin1=np.array([1+0j,0+1j],dtype=np.complex128)
        spin2=np.array([1+0j,1+0j],dtype=np.complex128)                
        spin1 /= np.linalg.norm(spin1)
        spin2 /= np.linalg.norm(spin2)
        print("initial spin 1",np.round(getBlochVector(spin1),3))
        print("initial spin 2",np.round(getBlochVector(spin2),3))
                        
        DT = float(self.txDT.text())
        print("initializinf dirac")
        
        self.dirac.initialize(L,DT,pos1,pos2,k1,k2,spin1,spin2)
        print("starting animation")
        self.canvasplot.im=None
        thmakepsi = Thread(target=self.threadPlayAnim,args=(self.dirac,))
        thmakepsi.start()
        
        
    def threadPlayAnim(self, objdirac):
        self.btPlay.setText("Stop")
        QtGui.QGuiApplication.processEvents()
        i=0
        while(self.busy):            
            print("steps",i)
            sex=objdirac.getSpinExpecValue()
            print("spin1",np.round(sex[0],4))
            print("spin2",np.round(sex[1],4))
            if i > 0:            
                objdirac.dostep()
            p= objdirac.getProb()
            self.canvasplot.plot(p)
            QtGui.QGuiApplication.processEvents()
            i+=1
            
            
        self.busy=False
        self.btPlay.setText("Play")    
        

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())

