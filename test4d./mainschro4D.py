import numpy as np
from time import time,sleep
from threading import Thread
from splitstep import SplitStepMethod
from PyQt5 import QtCore, QtGui, QtWidgets
from mainwnd import Ui_MainWindow
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from auxfunctions import make2DGaussian,makeGaussian2P2D,makeGaussian2P2Dold
from potentials import makeYBarrier
from schrod4D import schrod4D

class PlotCanvas(FigureCanvas):
    def __init__(self, parent = None, width = 5, height = 5, dpi = 100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(1,1,1)
 
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        self.figure.tight_layout()
        self.im=None
 
 
    def plot(self,p1,p2):            
        if self.im is None:
            self.im = self.axes.imshow(p1+p2,
                                   interpolation='bilinear',
                                   origin='lower', cmap='hot')
        else:                    
            self.im.set_data(p1+p2)
        self.draw()        
        
        
 
class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent=parent)
        self.setupUi(self)
        
        self.busy=False
        self.btPlay.clicked.connect(self.funcPlay)
        self.canvasplot = PlotCanvas(self.lbogl, width=4.5, height=4.5)
        self.canvasplot.move(0,0)
        self.Vplot=None
        self.schrod=None
        
    def funcPlay(self):
        self.busy=not self.busy        
        QtGui.QGuiApplication.processEvents()
        
        if not self.busy: return
        N=int(self.txN.text())
        L=float(self.txL.text())
        DT = float(self.txDT.text())
        pos1=np.array([float(self.txX1.text()),float(self.txY1.text())])
        pos2=np.array([float(self.txX2.text()),float(self.txY2.text())])
                       
        k1 = np.array([float(self.txKX1.text()),float(self.txKY1.text()),0.])
        k2 = np.array([float(self.txKX2.text()),float(self.txKY2.text()),0.])
        
        psi1 = make2DGaussian(N,L,k1,pos1)
        psi2 = make2DGaussian(N,L,k2,pos2)

        self.schrod=schrod4D(N)
        
        self.schrod.initialize(L,DT,pos1,pos2,k1,k2)
        #V = makeYBarrier(N,width=4,maxval=1e32)#np.zeros([N,N,N,N])
        #self.Vplot= np.einsum('jqik->ij', V)        
        #self.Vplot /=np.max(self.Vplot)        
                        
        
        # animation
        
        print("starting animation")
        self.canvasplot.im=None
        thmakepsi = Thread(target=self.threadPlayAnim,args=(self.schrod,))
        thmakepsi.start()
        
    def threadPlayAnim(self, schrod):
        self.btPlay.setText("Stop")
        QtGui.QGuiApplication.processEvents()
        i=0
        while(self.busy):            
            print("steps",i)            
            if i > 0:            
                schrod.dostep()
            
            p1,p2= schrod.getProb()
            # normalize so both plots have the same strenght()
            max1 = np.max(p1)
            max2 = np.max(p2)
            print("max psi 1",max1)
            print("max psi 2",max2)
            p1/=max1
            p2/=max2
            self.canvasplot.plot(p1,p2)
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

