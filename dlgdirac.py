# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dlgdirac.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_DlgDiract(object):
    def setupUi(self, DlgDiract):
        DlgDiract.setObjectName("DlgDiract")
        DlgDiract.resize(575, 546)
        font = QtGui.QFont()
        font.setPointSize(13)
        DlgDiract.setFont(font)
        self.buttonBox = QtWidgets.QDialogButtonBox(DlgDiract)
        self.buttonBox.setGeometry(QtCore.QRect(380, 440, 161, 81))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.buttonBox.setFont(font)
        self.buttonBox.setOrientation(QtCore.Qt.Vertical)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.groupBoxPos = QtWidgets.QGroupBox(DlgDiract)
        self.groupBoxPos.setGeometry(QtCore.QRect(20, 20, 191, 121))
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
        self.groupBox_K = QtWidgets.QGroupBox(DlgDiract)
        self.groupBox_K.setGeometry(QtCore.QRect(220, 20, 161, 131))
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
        self.groupBoxSpin = QtWidgets.QGroupBox(DlgDiract)
        self.groupBoxSpin.setGeometry(QtCore.QRect(0, 150, 581, 281))
        self.groupBoxSpin.setObjectName("groupBoxSpin")
        self.opcspin0 = QtWidgets.QRadioButton(self.groupBoxSpin)
        self.opcspin0.setGeometry(QtCore.QRect(10, 40, 311, 23))
        self.opcspin0.setChecked(False)
        self.opcspin0.setObjectName("opcspin0")
        self.opcspin1 = QtWidgets.QRadioButton(self.groupBoxSpin)
        self.opcspin1.setGeometry(QtCore.QRect(10, 70, 331, 23))
        self.opcspin1.setChecked(True)
        self.opcspin1.setObjectName("opcspin1")
        self.txspinupreal = QtWidgets.QLineEdit(self.groupBoxSpin)
        self.txspinupreal.setEnabled(True)
        self.txspinupreal.setGeometry(QtCore.QRect(100, 130, 121, 25))
        self.txspinupreal.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txspinupreal.setReadOnly(False)
        self.txspinupreal.setObjectName("txspinupreal")
        self.label_7 = QtWidgets.QLabel(self.groupBoxSpin)
        self.label_7.setGeometry(QtCore.QRect(130, 110, 71, 20))
        self.label_7.setObjectName("label_7")
        self.txspinupimag = QtWidgets.QLineEdit(self.groupBoxSpin)
        self.txspinupimag.setEnabled(True)
        self.txspinupimag.setGeometry(QtCore.QRect(230, 130, 121, 25))
        self.txspinupimag.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txspinupimag.setObjectName("txspinupimag")
        self.label_8 = QtWidgets.QLabel(self.groupBoxSpin)
        self.label_8.setGeometry(QtCore.QRect(250, 110, 111, 20))
        self.label_8.setObjectName("label_8")
        self.label_9 = QtWidgets.QLabel(self.groupBoxSpin)
        self.label_9.setGeometry(QtCore.QRect(20, 130, 71, 20))
        self.label_9.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_9.setObjectName("label_9")
        self.txspindownreal = QtWidgets.QLineEdit(self.groupBoxSpin)
        self.txspindownreal.setEnabled(True)
        self.txspindownreal.setGeometry(QtCore.QRect(100, 170, 121, 25))
        self.txspindownreal.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txspindownreal.setReadOnly(False)
        self.txspindownreal.setObjectName("txspindownreal")
        self.label_10 = QtWidgets.QLabel(self.groupBoxSpin)
        self.label_10.setGeometry(QtCore.QRect(0, 170, 91, 20))
        self.label_10.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_10.setObjectName("label_10")
        self.txspindownimag = QtWidgets.QLineEdit(self.groupBoxSpin)
        self.txspindownimag.setEnabled(True)
        self.txspindownimag.setGeometry(QtCore.QRect(230, 170, 121, 25))
        self.txspindownimag.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txspindownimag.setObjectName("txspindownimag")
        self.btPlotSpin = QtWidgets.QPushButton(self.groupBoxSpin)
        self.btPlotSpin.setEnabled(True)
        self.btPlotSpin.setGeometry(QtCore.QRect(120, 210, 181, 41))
        self.btPlotSpin.setObjectName("btPlotSpin")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.groupBoxSpin)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(380, 50, 191, 201))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.oglblochlayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.oglblochlayout.setContentsMargins(0, 0, 0, 0)
        self.oglblochlayout.setObjectName("oglblochlayout")
        self.groupBox_2 = QtWidgets.QGroupBox(DlgDiract)
        self.groupBox_2.setGeometry(QtCore.QRect(400, 20, 171, 131))
        font = QtGui.QFont()
        font.setFamily("Sans")
        font.setPointSize(13)
        font.setBold(True)
        font.setWeight(75)
        self.groupBox_2.setFont(font)
        self.groupBox_2.setStyleSheet("background-color: rgb(64, 98, 191);\n"
"color: rgb(238, 238, 236);")
        self.groupBox_2.setObjectName("groupBox_2")
        self.opcSplit = QtWidgets.QRadioButton(self.groupBox_2)
        self.opcSplit.setGeometry(QtCore.QRect(30, 40, 112, 23))
        font = QtGui.QFont()
        font.setPointSize(15)
        self.opcSplit.setFont(font)
        self.opcSplit.setObjectName("opcSplit")
        self.opcEigen = QtWidgets.QRadioButton(self.groupBox_2)
        self.opcEigen.setGeometry(QtCore.QRect(30, 80, 112, 23))
        font = QtGui.QFont()
        font.setPointSize(15)
        self.opcEigen.setFont(font)
        self.opcEigen.setChecked(True)
        self.opcEigen.setObjectName("opcEigen")
        self.label_12 = QtWidgets.QLabel(DlgDiract)
        self.label_12.setGeometry(QtCore.QRect(50, 450, 181, 31))
        self.label_12.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_12.setObjectName("label_12")
        self.txMaxSteps = QtWidgets.QLineEdit(DlgDiract)
        self.txMaxSteps.setGeometry(QtCore.QRect(250, 450, 61, 25))
        self.txMaxSteps.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txMaxSteps.setObjectName("txMaxSteps")

        self.retranslateUi(DlgDiract)
        self.buttonBox.accepted.connect(DlgDiract.accept)
        self.buttonBox.rejected.connect(DlgDiract.reject)
        self.opcspin0.toggled['bool'].connect(self.txspinupreal.setDisabled)
        self.opcspin0.toggled['bool'].connect(self.txspinupimag.setDisabled)
        self.opcspin0.toggled['bool'].connect(self.txspindownreal.setDisabled)
        self.opcspin0.toggled['bool'].connect(self.txspindownimag.setDisabled)
        self.opcspin0.toggled['bool'].connect(self.btPlotSpin.setDisabled)
        QtCore.QMetaObject.connectSlotsByName(DlgDiract)

    def retranslateUi(self, DlgDiract):
        _translate = QtCore.QCoreApplication.translate
        DlgDiract.setWindowTitle(_translate("DlgDiract", "Dirac options"))
        self.groupBoxPos.setTitle(_translate("DlgDiract", "inital gaussian position"))
        self.label_3.setText(_translate("DlgDiract", "Z"))
        self.label_2.setText(_translate("DlgDiract", "Y"))
        self.label.setText(_translate("DlgDiract", "X"))
        self.txposX.setToolTip(_translate("DlgDiract", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txposX.setText(_translate("DlgDiract", "-0.7"))
        self.txposY.setToolTip(_translate("DlgDiract", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txposY.setText(_translate("DlgDiract", "0.0"))
        self.txposZ.setToolTip(_translate("DlgDiract", "<html><head/><body><p>from -0.8 to 0.8</p></body></html>"))
        self.txposZ.setText(_translate("DlgDiract", "0.0"))
        self.groupBox_K.setToolTip(_translate("DlgDiract", "<html><head/><body><p>spatial frequency with respect to spatial extent of the simulation</p><p>p = k * np.pi / L</p></body></html>"))
        self.groupBox_K.setTitle(_translate("DlgDiract", "K  (gaussian speed)"))
        self.label_4.setText(_translate("DlgDiract", "k Z"))
        self.label_5.setText(_translate("DlgDiract", "k Y"))
        self.label_6.setText(_translate("DlgDiract", "k X"))
        self.txkX.setText(_translate("DlgDiract", "22"))
        self.txkY.setText(_translate("DlgDiract", "0"))
        self.txkZ.setText(_translate("DlgDiract", "0"))
        self.groupBoxSpin.setTitle(_translate("DlgDiract", "inital Spin and spinors"))
        self.opcspin0.setToolTip(_translate("DlgDiract", "<html><head/><body><p>The spinors will be created from momentun with spin up in Z axis</p></body></html>"))
        self.opcspin0.setText(_translate("DlgDiract", "spin and spinors from kinetic"))
        self.opcspin1.setToolTip(_translate("DlgDiract", "<html><head/><body><p>specify the spin you want below amd the spinors will be created from this and energy eigen Spinors mixing your spin and the momentum</p></body></html>"))
        self.opcspin1.setText(_translate("DlgDiract", "make spinors from specific spin below:"))
        self.txspinupreal.setText(_translate("DlgDiract", "0.5"))
        self.label_7.setText(_translate("DlgDiract", "real"))
        self.txspinupimag.setText(_translate("DlgDiract", "0.3"))
        self.label_8.setText(_translate("DlgDiract", "imaginary"))
        self.label_9.setText(_translate("DlgDiract", "spin up"))
        self.txspindownreal.setText(_translate("DlgDiract", "0.0"))
        self.label_10.setText(_translate("DlgDiract", "spin down"))
        self.txspindownimag.setText(_translate("DlgDiract", "0.5"))
        self.btPlotSpin.setToolTip(_translate("DlgDiract", "<html><head/><body><p>renders the spin normalized with OpenGL in the small window at the right</p></body></html>"))
        self.btPlotSpin.setText(_translate("DlgDiract", "Plot spin"))
        self.groupBox_2.setTitle(_translate("DlgDiract", "Simulation Mode "))
        self.opcSplit.setText(_translate("DlgDiract", "SplitStep"))
        self.opcEigen.setText(_translate("DlgDiract", "free eigen"))
        self.label_12.setText(_translate("DlgDiract", "max animation frames "))
        self.txMaxSteps.setToolTip(_translate("DlgDiract", "<html><head/><body><p>max numer of frames animation ( you can stop before )</p></body></html>"))
        self.txMaxSteps.setText(_translate("DlgDiract", "150"))

