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
        MainWindow.resize(556, 678)
        font = QtGui.QFont()
        font.setPointSize(13)
        MainWindow.setFont(font)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.groupBoxPos = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBoxPos.setGeometry(QtCore.QRect(50, 10, 121, 151))
        self.groupBoxPos.setObjectName("groupBoxPos")
        self.label_2 = QtWidgets.QLabel(self.groupBoxPos)
        self.label_2.setGeometry(QtCore.QRect(20, 60, 21, 17))
        self.label_2.setObjectName("label_2")
        self.label = QtWidgets.QLabel(self.groupBoxPos)
        self.label.setGeometry(QtCore.QRect(20, 30, 21, 17))
        self.label.setObjectName("label")
        self.txX1 = QtWidgets.QLineEdit(self.groupBoxPos)
        self.txX1.setGeometry(QtCore.QRect(40, 30, 51, 25))
        self.txX1.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txX1.setObjectName("txX1")
        self.txY1 = QtWidgets.QLineEdit(self.groupBoxPos)
        self.txY1.setGeometry(QtCore.QRect(40, 60, 51, 25))
        self.txY1.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txY1.setObjectName("txY1")
        self.label_5 = QtWidgets.QLabel(self.groupBoxPos)
        self.label_5.setGeometry(QtCore.QRect(10, 120, 31, 20))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.groupBoxPos)
        self.label_6.setGeometry(QtCore.QRect(10, 90, 31, 20))
        self.label_6.setObjectName("label_6")
        self.txKY1 = QtWidgets.QLineEdit(self.groupBoxPos)
        self.txKY1.setGeometry(QtCore.QRect(40, 120, 61, 25))
        self.txKY1.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txKY1.setObjectName("txKY1")
        self.txKX1 = QtWidgets.QLineEdit(self.groupBoxPos)
        self.txKX1.setGeometry(QtCore.QRect(40, 90, 61, 25))
        self.txKX1.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txKX1.setObjectName("txKX1")
        self.btPlay = QtWidgets.QPushButton(self.centralwidget)
        self.btPlay.setGeometry(QtCore.QRect(450, 60, 89, 51))
        font = QtGui.QFont()
        font.setPointSize(15)
        self.btPlay.setFont(font)
        self.btPlay.setObjectName("btPlay")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(340, 40, 31, 20))
        self.label_7.setObjectName("label_7")
        self.txN = QtWidgets.QLineEdit(self.centralwidget)
        self.txN.setGeometry(QtCore.QRect(360, 40, 41, 25))
        self.txN.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txN.setObjectName("txN")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(340, 80, 31, 20))
        self.label_8.setObjectName("label_8")
        self.txL = QtWidgets.QLineEdit(self.centralwidget)
        self.txL.setGeometry(QtCore.QRect(360, 80, 41, 25))
        self.txL.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txL.setObjectName("txL")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(330, 120, 31, 20))
        self.label_9.setObjectName("label_9")
        self.txDT = QtWidgets.QLineEdit(self.centralwidget)
        self.txDT.setGeometry(QtCore.QRect(360, 120, 81, 25))
        self.txDT.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txDT.setObjectName("txDT")
        self.groupBoxPos_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBoxPos_2.setGeometry(QtCore.QRect(190, 10, 121, 151))
        self.groupBoxPos_2.setObjectName("groupBoxPos_2")
        self.label_3 = QtWidgets.QLabel(self.groupBoxPos_2)
        self.label_3.setGeometry(QtCore.QRect(20, 60, 21, 17))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.groupBoxPos_2)
        self.label_4.setGeometry(QtCore.QRect(20, 30, 21, 17))
        self.label_4.setObjectName("label_4")
        self.txX2 = QtWidgets.QLineEdit(self.groupBoxPos_2)
        self.txX2.setGeometry(QtCore.QRect(40, 30, 51, 25))
        self.txX2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txX2.setObjectName("txX2")
        self.txY2 = QtWidgets.QLineEdit(self.groupBoxPos_2)
        self.txY2.setGeometry(QtCore.QRect(40, 60, 51, 25))
        self.txY2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txY2.setObjectName("txY2")
        self.label_10 = QtWidgets.QLabel(self.groupBoxPos_2)
        self.label_10.setGeometry(QtCore.QRect(10, 120, 31, 20))
        self.label_10.setObjectName("label_10")
        self.label_11 = QtWidgets.QLabel(self.groupBoxPos_2)
        self.label_11.setGeometry(QtCore.QRect(10, 90, 31, 20))
        self.label_11.setObjectName("label_11")
        self.txKY2 = QtWidgets.QLineEdit(self.groupBoxPos_2)
        self.txKY2.setGeometry(QtCore.QRect(40, 120, 61, 25))
        self.txKY2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txKY2.setObjectName("txKY2")
        self.txKX2 = QtWidgets.QLineEdit(self.groupBoxPos_2)
        self.txKX2.setGeometry(QtCore.QRect(40, 90, 61, 25))
        self.txKX2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txKX2.setObjectName("txKX2")
        self.lbogl = QtWidgets.QLabel(self.centralwidget)
        self.lbogl.setGeometry(QtCore.QRect(20, 190, 451, 441))
        self.lbogl.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lbogl.setText("")
        self.lbogl.setObjectName("lbogl")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 556, 25))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.groupBoxPos.setTitle(_translate("MainWindow", "particle 1"))
        self.label_2.setText(_translate("MainWindow", "Y"))
        self.label.setText(_translate("MainWindow", "X"))
        self.txX1.setToolTip(_translate("MainWindow", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txX1.setText(_translate("MainWindow", "-0.5"))
        self.txY1.setToolTip(_translate("MainWindow", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txY1.setText(_translate("MainWindow", "0"))
        self.label_5.setText(_translate("MainWindow", "k Y"))
        self.label_6.setText(_translate("MainWindow", "k X"))
        self.txKY1.setText(_translate("MainWindow", "0"))
        self.txKX1.setText(_translate("MainWindow", "-10"))
        self.btPlay.setText(_translate("MainWindow", "Play"))
        self.label_7.setText(_translate("MainWindow", "N"))
        self.txN.setText(_translate("MainWindow", "32"))
        self.label_8.setText(_translate("MainWindow", "L"))
        self.txL.setText(_translate("MainWindow", "200"))
        self.label_9.setText(_translate("MainWindow", "DT"))
        self.txDT.setText(_translate("MainWindow", "20"))
        self.groupBoxPos_2.setTitle(_translate("MainWindow", "particle 2"))
        self.label_3.setText(_translate("MainWindow", "Y"))
        self.label_4.setText(_translate("MainWindow", "X"))
        self.txX2.setToolTip(_translate("MainWindow", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txX2.setText(_translate("MainWindow", "0"))
        self.txY2.setToolTip(_translate("MainWindow", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txY2.setText(_translate("MainWindow", "-0.5"))
        self.label_10.setText(_translate("MainWindow", "k Y"))
        self.label_11.setText(_translate("MainWindow", "k X"))
        self.txKY2.setText(_translate("MainWindow", "10"))
        self.txKX2.setText(_translate("MainWindow", "0"))

