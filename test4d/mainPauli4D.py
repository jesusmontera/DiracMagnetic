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
        self.btLoadSpin.clicked.connect(self.funcLoadSpin)
        self.btReplay.clicked.connect(self.funcReplay)
        self.canvasplot = PlotCanvas(self.lbogl, width=4.5, height=4.5)
        self.canvasplot.move(0,0)
        self.pauli=None                
        self.framesanim=None
        spin=np.load("spins_Yrightright_noentangled.npy")
        self.spintotextboxes(spin)
                
                
    def spintotextboxes(self,spin):
        self.txspin0r.setText('{:.6e}'.format(spin[0].real))
        self.txspin0i.setText('{:.6e}'.format(spin[0].imag))
        self.txspin1r.setText('{:.6e}'.format(spin[1].real))
        self.txspin1i.setText('{:.6e}'.format(spin[1].imag))
        self.txspin2r.setText('{:.6e}'.format(spin[2].real))
        self.txspin2i.setText('{:.6e}'.format(spin[2].imag))
        self.txspin3r.setText('{:.6e}'.format(spin[3].real))
        self.txspin3i.setText('{:.6e}'.format(spin[3].imag))
        self.spinparse()
        
    def funcLoadSpin(self):
        file , check = QtWidgets.QFileDialog.getOpenFileName(self,                                                             
                                                             "Open spin npy file",
                                                             "", "npy Files (*.npy)")        
        
        if check:                        
            spin= np.load(file)
            self.spintotextboxes(spin)            
            
    def spinparse(self):
        
        c0=complex(float(self.txspin0r.text()) + float(self.txspin0i.text())*1j)
        c1=complex(float(self.txspin1r.text()) + float(self.txspin1i.text())*1j)
        c2=complex(float(self.txspin2r.text()) + float(self.txspin2i.text())*1j)
        c3=complex(float(self.txspin3r.text()) + float(self.txspin3i.text())*1j)
        
        spin= np.array([c0,c1,c2,c3])
        spin /= np.linalg.norm(spin)
        
        bloch1=bloch_vector_from_state_vector(spin,0)
        bloch2=bloch_vector_from_state_vector(spin,1)
        # spin entaglement quantification
        rho1= density_matrix_from_state_vector(spin,[0])
        rho2=  density_matrix_from_state_vector(spin,[1])
        entropy1 = np.round(von_newman_entropy(rho1),4)
        entropy2 = np.round(von_newman_entropy(rho2),4)        

        self.lbbloch1.setText("bloch 1:  X " + str(np.round(bloch1[0],5))+
                              "  Y " + str(np.round(bloch1[1],5))+
                              "  Z " + str(np.round(bloch1[2],5)) + " entagled="+str(entropy1))
        self.lbbloch2.setText("bloch 2:  X " + str(np.round(bloch2[0],5))+
                              "  Y " + str(np.round(bloch2[1],5))+
                              "  Z " + str(np.round(bloch2[2],5))+ " entagled="+str(entropy2))
        return spin
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
        
        if not self.busy: return
        self.btPlay.setText("Stop")
        QtGui.QGuiApplication.processEvents()
        N=int(self.txN.text())
        L=float(self.txL.text())                
        self.pauli= Pauli4D()        
                
        pos=np.array([float(self.txX1.text()),float(self.txY1.text()),
                      float(self.txX2.text()),float(self.txY2.text())])
        
                       
        k=np.array([float(self.txKX1.text()),float(self.txKY1.text()),                    
                    float(self.txKX2.text()),float(self.txKY2.text())])
        
        # load the initial spin state vector(4 complex numbers) that can be entangled or not
        # remove this and specify the one you like
        


                
        DT = float(self.txDT.text())
        print("initializing Pauli...")
        spin = self.spinparse()        
        self.pauli.initialize(N,L,DT,pos,k,spin)
        self.pauli.makeBspinnors(self.radioBhalf.isChecked())# make B spinnor potential matrix

        print("starting animation")
        self.canvasplot.im=None
        thmakepsi = Thread(target=self.threadPlayAnim,args=(self.pauli,))
        thmakepsi.start()
        
        
    def threadPlayAnim(self, obj):        
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

