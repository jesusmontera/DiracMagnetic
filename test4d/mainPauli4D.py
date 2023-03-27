import numpy as np
from time import time,sleep
from threading import Thread
from PyQt5 import QtCore, QtGui, QtWidgets
from mainwnd import Ui_MainWindow
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from Pauli4D import Pauli4D
from auxfunctions import bloch_vector_from_state_vector,von_newman_entropy,density_matrix_from_state_vector

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
        self.btReplay.clicked.connect(self.funcReplay)
        self.canvasplot = PlotCanvas(self.lbogl, width=4.5, height=4.5)
        self.canvasplot.move(0,0)
        self.pauli=None                
        self.framesanim=None
    def funcReplay(self):
        if self.framesanim is None: return
        self.busy=not self.busy
        if self.busy:           
            n= len(self.framesanim)
            self.btReplay.setText("stop")
            for i in range(n):
                self.canvasplot.plot(self.framesanim[i])
                tini=time()
                while self.busy:
                    QtGui.QGuiApplication.processEvents()
                    if time()-tini>0.1:break
                if not self.busy:
                    break
        self.btReplay.setText("Replay")
        self.busy=False
        
            
        print("hola")
    def funcPlay(self):
        self.busy=not self.busy        
        QtGui.QGuiApplication.processEvents()
        
        if not self.busy: return
        
        N=int(self.txN.text())
        L=float(self.txL.text())                
        self.pauli= Pauli4D()        
                
        pos=np.array([float(self.txX1.text()),float(self.txY1.text()),
                      float(self.txX2.text()),float(self.txY2.text())])
        
                       
        k=np.array([float(self.txKX1.text()),float(self.txKY1.text()),                    
                    float(self.txKX2.text()),float(self.txKY2.text())])
        
        # load the initial spin state vector(4 complex numbers) that can be entangled or not
        # remove this and specify the one you like
        initial_spin=np.load("spinsentangled.npy") 
        print("initial spin ",np.round(initial_spin,3))                
        bloch0=bloch_vector_from_state_vector(initial_spin,0)
        bloch1=bloch_vector_from_state_vector(initial_spin,1)
        print("initial spin bloch 1",np.round(bloch0,3))
        print("initial spin bloch 2",np.round(bloch1,3))        
        
        # print spin entaglement quantification
        rho1= density_matrix_from_state_vector(initial_spin,[0])
        rho2=  density_matrix_from_state_vector(initial_spin,[1])
        print("initial spin 1 entropy (entangled) = ",np.round(von_newman_entropy(rho1.data),6))
        print("initial spin 2 entropy (entangled) = ",np.round(von_newman_entropy(rho2.data),6))

                
        DT = float(self.txDT.text())
        print("initializing Pauli...")
        
        self.pauli.initialize(N,L,DT,pos,k,initial_spin)
        

        print("starting animation")
        self.canvasplot.im=None
        thmakepsi = Thread(target=self.threadPlayAnim,args=(self.pauli,))
        thmakepsi.start()
        
        
    def threadPlayAnim(self, obj):
        self.btPlay.setText("Stop")
        QtGui.QGuiApplication.processEvents()
        i=0
        self.framesanim=[]
        while(self.busy):            
            print("steps",i)            

            if i > 0:            
                obj.dostep()

            sex=self.pauli.getSpinExpecValue()            
            print("spin1",np.round(sex[0],4))
            print("spin2",np.round(sex[1],4))
            
            p= obj.getProb()
            self.framesanim.append(p)
            self.canvasplot.plot(p)
            QtGui.QGuiApplication.processEvents()
            i+=1
            if i>40: break
                
            
            
        self.busy=False
        self.btPlay.setText("Play")    
        

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())

