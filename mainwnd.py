# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainwnd.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(725, 673)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.btPlay = QtWidgets.QPushButton(self.centralwidget)
        self.btPlay.setEnabled(False)
        self.btPlay.setGeometry(QtCore.QRect(520, 20, 131, 61))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.btPlay.setFont(font)
        self.btPlay.setObjectName("btPlay")
        self.lbinfo = QtWidgets.QLabel(self.centralwidget)
        self.lbinfo.setGeometry(QtCore.QRect(40, 100, 571, 41))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.lbinfo.setFont(font)
        self.lbinfo.setObjectName("lbinfo")
        self.txN = QtWidgets.QLineEdit(self.centralwidget)
        self.txN.setEnabled(True)
        self.txN.setGeometry(QtCore.QRect(330, 10, 51, 25))
        font = QtGui.QFont()
        font.setPointSize(13)
        self.txN.setFont(font)
        self.txN.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txN.setReadOnly(True)
        self.txN.setObjectName("txN")
        self.sliderOGL = QtWidgets.QSlider(self.centralwidget)
        self.sliderOGL.setGeometry(QtCore.QRect(10, 130, 641, 31))
        font = QtGui.QFont()
        font.setPointSize(13)
        self.sliderOGL.setFont(font)
        self.sliderOGL.setOrientation(QtCore.Qt.Horizontal)
        self.sliderOGL.setObjectName("sliderOGL")
        self.btMake = QtWidgets.QPushButton(self.centralwidget)
        self.btMake.setGeometry(QtCore.QRect(20, 20, 131, 61))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.btMake.setFont(font)
        self.btMake.setObjectName("btMake")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(310, 10, 21, 17))
        self.label.setObjectName("label")
        self.txL = QtWidgets.QLineEdit(self.centralwidget)
        self.txL.setEnabled(True)
        self.txL.setGeometry(QtCore.QRect(330, 44, 51, 21))
        font = QtGui.QFont()
        font.setPointSize(13)
        self.txL.setFont(font)
        self.txL.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txL.setReadOnly(False)
        self.txL.setObjectName("txL")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(310, 40, 21, 17))
        self.label_2.setObjectName("label_2")
        self.txDT = QtWidgets.QLineEdit(self.centralwidget)
        self.txDT.setGeometry(QtCore.QRect(400, 30, 81, 25))
        font = QtGui.QFont()
        font.setPointSize(13)
        self.txDT.setFont(font)
        self.txDT.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txDT.setObjectName("txDT")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(410, 10, 31, 17))
        self.label_3.setObjectName("label_3")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(10, 170, 661, 481))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.oglalyout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.oglalyout.setContentsMargins(0, 0, 0, 0)
        self.oglalyout.setObjectName("oglalyout")
        self.radioDirac = QtWidgets.QRadioButton(self.centralwidget)
        self.radioDirac.setGeometry(QtCore.QRect(170, 50, 71, 23))
        font = QtGui.QFont()
        font.setPointSize(13)
        self.radioDirac.setFont(font)
        self.radioDirac.setChecked(False)
        self.radioDirac.setObjectName("radioDirac")
        self.radioSchro = QtWidgets.QRadioButton(self.centralwidget)
        self.radioSchro.setGeometry(QtCore.QRect(170, 80, 141, 23))
        font = QtGui.QFont()
        font.setPointSize(13)
        self.radioSchro.setFont(font)
        self.radioSchro.setObjectName("radioSchro")
        self.radioPauli = QtWidgets.QRadioButton(self.centralwidget)
        self.radioPauli.setGeometry(QtCore.QRect(170, 20, 141, 23))
        font = QtGui.QFont()
        font.setPointSize(13)
        self.radioPauli.setFont(font)
        self.radioPauli.setChecked(True)
        self.radioPauli.setObjectName("radioPauli")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 725, 22))
        self.menubar.setObjectName("menubar")
        self.menumi_menu = QtWidgets.QMenu(self.menubar)
        self.menumi_menu.setObjectName("menumi_menu")
        self.menuLoad = QtWidgets.QMenu(self.menumi_menu)
        self.menuLoad.setObjectName("menuLoad")
        self.menuSave = QtWidgets.QMenu(self.menumi_menu)
        self.menuSave.setObjectName("menuSave")
        self.menuMagnetic_Field = QtWidgets.QMenu(self.menubar)
        self.menuMagnetic_Field.setObjectName("menuMagnetic_Field")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.mnuLoadB = QtWidgets.QAction(MainWindow)
        self.mnuLoadB.setObjectName("mnuLoadB")
        self.mnuLoadSch = QtWidgets.QAction(MainWindow)
        self.mnuLoadSch.setObjectName("mnuLoadSch")
        self.mnusaveB = QtWidgets.QAction(MainWindow)
        self.mnusaveB.setObjectName("mnusaveB")
        self.mnuSaveSch = QtWidgets.QAction(MainWindow)
        self.mnuSaveSch.setObjectName("mnuSaveSch")
        self.mnushowSch = QtWidgets.QAction(MainWindow)
        self.mnushowSch.setObjectName("mnushowSch")
        self.mnuSaveDirac = QtWidgets.QAction(MainWindow)
        self.mnuSaveDirac.setObjectName("mnuSaveDirac")
        self.mnuLoadDirac = QtWidgets.QAction(MainWindow)
        self.mnuLoadDirac.setObjectName("mnuLoadDirac")
        self.mnumagshow = QtWidgets.QAction(MainWindow)
        self.mnumagshow.setObjectName("mnumagshow")
        self.menuLoad.addAction(self.mnuLoadSch)
        self.menuLoad.addAction(self.mnuLoadDirac)
        self.menuSave.addAction(self.mnuSaveSch)
        self.menuSave.addAction(self.mnuSaveDirac)
        self.menumi_menu.addAction(self.menuLoad.menuAction())
        self.menumi_menu.addAction(self.menuSave.menuAction())
        self.menuMagnetic_Field.addAction(self.mnumagshow)
        self.menubar.addAction(self.menumi_menu.menuAction())
        self.menubar.addAction(self.menuMagnetic_Field.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Dirac magnetic field in OpenGL"))
        self.btPlay.setText(_translate("MainWindow", "Play OpenGL"))
        self.lbinfo.setText(_translate("MainWindow", "info"))
        self.txN.setToolTip(_translate("MainWindow", "<html><head/><body><p>Number of points for PSI 3d </p><p>array[N,N,N]</p></body></html>"))
        self.txN.setText(_translate("MainWindow", "64"))
        self.sliderOGL.setToolTip(_translate("MainWindow", "animation position "))
        self.btMake.setToolTip(_translate("MainWindow", "Select Dirac or Scroduinger  and make the psi animation"))
        self.btMake.setText(_translate("MainWindow", "Make PSI"))
        self.label.setText(_translate("MainWindow", "N"))
        self.txL.setToolTip(_translate("MainWindow", "<html><head/><body><p>Box space size in AU(Hartree atomic units=Bohrn length)</p></body></html>"))
        self.txL.setText(_translate("MainWindow", "2.0"))
        self.label_2.setText(_translate("MainWindow", "L"))
        self.txDT.setToolTip(_translate("MainWindow", "<html><head/><body><p>time step.</p><p>if you decrease kinetic you should decrease DT</p></body></html>"))
        self.txDT.setText(_translate("MainWindow", "0.001"))
        self.label_3.setText(_translate("MainWindow", "DT"))
        self.radioDirac.setText(_translate("MainWindow", "Dirac"))
        self.radioSchro.setText(_translate("MainWindow", "Schroduinger"))
        self.radioPauli.setText(_translate("MainWindow", "Pauli"))
        self.menumi_menu.setTitle(_translate("MainWindow", "File"))
        self.menuLoad.setTitle(_translate("MainWindow", "Load"))
        self.menuSave.setTitle(_translate("MainWindow", "Save"))
        self.menuMagnetic_Field.setTitle(_translate("MainWindow", "Magnetic Field"))
        self.mnuLoadB.setText(_translate("MainWindow", "magnetic field"))
        self.mnuLoadSch.setText(_translate("MainWindow", "schroduinger"))
        self.mnusaveB.setText(_translate("MainWindow", "magnetic field"))
        self.mnuSaveSch.setText(_translate("MainWindow", "Schroduinger"))
        self.mnushowSch.setText(_translate("MainWindow", "Show window"))
        self.mnuSaveDirac.setText(_translate("MainWindow", "Dirac"))
        self.mnuLoadDirac.setText(_translate("MainWindow", "Dirac"))
        self.mnumagshow.setText(_translate("MainWindow", "Show window"))

