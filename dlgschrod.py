# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dlgschrod.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_DlgSch(object):
    def setupUi(self, DlgSch):
        DlgSch.setObjectName("DlgSch")
        DlgSch.resize(400, 418)
        font = QtGui.QFont()
        font.setPointSize(13)
        DlgSch.setFont(font)
        self.buttonBox = QtWidgets.QDialogButtonBox(DlgSch)
        self.buttonBox.setGeometry(QtCore.QRect(30, 350, 241, 61))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.groupBoxPos = QtWidgets.QGroupBox(DlgSch)
        self.groupBoxPos.setGeometry(QtCore.QRect(20, 20, 191, 131))
        self.groupBoxPos.setObjectName("groupBoxPos")
        self.label_3 = QtWidgets.QLabel(self.groupBoxPos)
        self.label_3.setGeometry(QtCore.QRect(30, 90, 21, 17))
        self.label_3.setObjectName("label_3")
        self.label_2 = QtWidgets.QLabel(self.groupBoxPos)
        self.label_2.setGeometry(QtCore.QRect(30, 60, 21, 17))
        self.label_2.setObjectName("label_2")
        self.label = QtWidgets.QLabel(self.groupBoxPos)
        self.label.setGeometry(QtCore.QRect(30, 30, 21, 17))
        self.label.setObjectName("label")
        self.txposX = QtWidgets.QLineEdit(self.groupBoxPos)
        self.txposX.setGeometry(QtCore.QRect(60, 30, 61, 25))
        self.txposX.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txposX.setObjectName("txposX")
        self.txposY = QtWidgets.QLineEdit(self.groupBoxPos)
        self.txposY.setGeometry(QtCore.QRect(60, 60, 61, 25))
        self.txposY.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txposY.setObjectName("txposY")
        self.txposZ = QtWidgets.QLineEdit(self.groupBoxPos)
        self.txposZ.setGeometry(QtCore.QRect(60, 90, 61, 25))
        self.txposZ.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txposZ.setObjectName("txposZ")
        self.groupBox_K = QtWidgets.QGroupBox(DlgSch)
        self.groupBox_K.setGeometry(QtCore.QRect(220, 20, 181, 131))
        self.groupBox_K.setObjectName("groupBox_K")
        self.label_4 = QtWidgets.QLabel(self.groupBox_K)
        self.label_4.setGeometry(QtCore.QRect(30, 90, 21, 17))
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.groupBox_K)
        self.label_5.setGeometry(QtCore.QRect(30, 60, 31, 20))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.groupBox_K)
        self.label_6.setGeometry(QtCore.QRect(30, 30, 31, 20))
        self.label_6.setObjectName("label_6")
        self.txkX = QtWidgets.QLineEdit(self.groupBox_K)
        self.txkX.setGeometry(QtCore.QRect(60, 30, 61, 25))
        self.txkX.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txkX.setObjectName("txkX")
        self.txkY = QtWidgets.QLineEdit(self.groupBox_K)
        self.txkY.setGeometry(QtCore.QRect(60, 60, 61, 25))
        self.txkY.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txkY.setObjectName("txkY")
        self.txkZ = QtWidgets.QLineEdit(self.groupBox_K)
        self.txkZ.setGeometry(QtCore.QRect(60, 90, 61, 25))
        self.txkZ.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txkZ.setObjectName("txkZ")
        self.groupBoxSpin = QtWidgets.QGroupBox(DlgSch)
        self.groupBoxSpin.setGeometry(QtCore.QRect(50, 160, 271, 91))
        self.groupBoxSpin.setObjectName("groupBoxSpin")
        self.label_7 = QtWidgets.QLabel(self.groupBoxSpin)
        self.label_7.setGeometry(QtCore.QRect(40, 30, 21, 17))
        self.label_7.setObjectName("label_7")
        self.txSpinX = QtWidgets.QLineEdit(self.groupBoxSpin)
        self.txSpinX.setGeometry(QtCore.QRect(20, 50, 61, 25))
        self.txSpinX.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txSpinX.setObjectName("txSpinX")
        self.label_8 = QtWidgets.QLabel(self.groupBoxSpin)
        self.label_8.setGeometry(QtCore.QRect(120, 30, 21, 17))
        self.label_8.setObjectName("label_8")
        self.txSpinY = QtWidgets.QLineEdit(self.groupBoxSpin)
        self.txSpinY.setGeometry(QtCore.QRect(100, 50, 61, 25))
        self.txSpinY.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txSpinY.setObjectName("txSpinY")
        self.txSpinZ = QtWidgets.QLineEdit(self.groupBoxSpin)
        self.txSpinZ.setGeometry(QtCore.QRect(170, 50, 61, 25))
        self.txSpinZ.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txSpinZ.setObjectName("txSpinZ")
        self.label_9 = QtWidgets.QLabel(self.groupBoxSpin)
        self.label_9.setGeometry(QtCore.QRect(190, 30, 21, 17))
        self.label_9.setObjectName("label_9")
        self.checkB = QtWidgets.QCheckBox(DlgSch)
        self.checkB.setGeometry(QtCore.QRect(30, 260, 201, 41))
        font = QtGui.QFont()
        font.setPointSize(13)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.checkB.setFont(font)
        self.checkB.setObjectName("checkB")
        self.label_11 = QtWidgets.QLabel(DlgSch)
        self.label_11.setGeometry(QtCore.QRect(230, 270, 151, 17))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.label_12 = QtWidgets.QLabel(DlgSch)
        self.label_12.setGeometry(QtCore.QRect(20, 310, 181, 31))
        self.label_12.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_12.setObjectName("label_12")
        self.txMaxSteps = QtWidgets.QLineEdit(DlgSch)
        self.txMaxSteps.setGeometry(QtCore.QRect(200, 310, 61, 25))
        self.txMaxSteps.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txMaxSteps.setObjectName("txMaxSteps")

        self.retranslateUi(DlgSch)
        self.buttonBox.accepted.connect(DlgSch.accept)
        self.buttonBox.rejected.connect(DlgSch.reject)
        QtCore.QMetaObject.connectSlotsByName(DlgSch)

    def retranslateUi(self, DlgSch):
        _translate = QtCore.QCoreApplication.translate
        DlgSch.setWindowTitle(_translate("DlgSch", "Schroduinger  options"))
        self.groupBoxPos.setTitle(_translate("DlgSch", "inital gaussian position"))
        self.label_3.setText(_translate("DlgSch", "Z"))
        self.label_2.setText(_translate("DlgSch", "Y"))
        self.label.setText(_translate("DlgSch", "X"))
        self.txposX.setToolTip(_translate("DlgSch", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txposX.setText(_translate("DlgSch", "-0.7"))
        self.txposY.setToolTip(_translate("DlgSch", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txposY.setText(_translate("DlgSch", "0.0"))
        self.txposZ.setToolTip(_translate("DlgSch", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txposZ.setText(_translate("DlgSch", "0.0"))
        self.groupBox_K.setToolTip(_translate("DlgSch", "<html><head/><body><p>momentum (speed) is:</p><p>p=k*2*pi / L </p></body></html>"))
        self.groupBox_K.setTitle(_translate("DlgSch", "kinetic Energy"))
        self.label_4.setText(_translate("DlgSch", "k Z"))
        self.label_5.setText(_translate("DlgSch", "k Y"))
        self.label_6.setText(_translate("DlgSch", "k X"))
        self.txkX.setText(_translate("DlgSch", "22"))
        self.txkY.setText(_translate("DlgSch", "0"))
        self.txkZ.setText(_translate("DlgSch", "0"))
        self.groupBoxSpin.setToolTip(_translate("DlgSch", "<html><head/><body><p>dir vector for spin that will be used to interact with magnetic field. WIll be automatically normalized</p></body></html>"))
        self.groupBoxSpin.setTitle(_translate("DlgSch", "initial Spin in Bloch coordinates"))
        self.label_7.setText(_translate("DlgSch", "X"))
        self.txSpinX.setToolTip(_translate("DlgSch", "<html><head/><body><p>X dir vector </p><p>from -1 to 1</p></body></html>"))
        self.txSpinX.setText(_translate("DlgSch", "1"))
        self.label_8.setText(_translate("DlgSch", "Y"))
        self.txSpinY.setToolTip(_translate("DlgSch", "<html><head/><body><p>X dir vector </p><p>from -1 to 1</p></body></html>"))
        self.txSpinY.setText(_translate("DlgSch", "0"))
        self.txSpinZ.setToolTip(_translate("DlgSch", "<html><head/><body><p>X dir vector </p><p>from -1 to 1</p></body></html>"))
        self.txSpinZ.setText(_translate("DlgSch", "0"))
        self.label_9.setText(_translate("DlgSch", "Z"))
        self.checkB.setToolTip(_translate("DlgSch", "<html><head/><body><p>the field created or loaded in the magnetic window will be appplied by dot spin and B after U step of DIrac</p></body></html>"))
        self.checkB.setText(_translate("DlgSch", "Apply magnetic B field"))
        self.label_11.setText(_translate("DlgSch", "(from magnetic dialog)"))
        self.label_12.setText(_translate("DlgSch", "max animation frames "))
        self.txMaxSteps.setToolTip(_translate("DlgSch", "<html><head/><body><p>max numer of frames animation ( you can stop before )</p></body></html>"))
        self.txMaxSteps.setText(_translate("DlgSch", "150"))

