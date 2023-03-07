# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dlgB.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(642, 653)
        font = QtGui.QFont()
        font.setPointSize(13)
        Dialog.setFont(font)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(510, 30, 121, 91))
        self.buttonBox.setOrientation(QtCore.Qt.Vertical)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(190, 10, 51, 20))
        self.label.setObjectName("label")
        self.txmaxB = QtWidgets.QLineEdit(Dialog)
        self.txmaxB.setGeometry(QtCore.QRect(170, 40, 81, 25))
        self.txmaxB.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.txmaxB.setObjectName("txmaxB")
        self.btMakeB = QtWidgets.QPushButton(Dialog)
        self.btMakeB.setGeometry(QtCore.QRect(390, 40, 89, 61))
        self.btMakeB.setObjectName("btMakeB")
        self.verticalLayoutWidget = QtWidgets.QWidget(Dialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(40, 160, 591, 481))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.oglayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.oglayout.setContentsMargins(0, 0, 0, 0)
        self.oglayout.setObjectName("oglayout")
        self.listdir = QtWidgets.QListWidget(Dialog)
        self.listdir.setGeometry(QtCore.QRect(270, 40, 121, 61))
        self.listdir.setObjectName("listdir")
        item = QtWidgets.QListWidgetItem()
        self.listdir.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.listdir.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.listdir.addItem(item)
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(270, 10, 161, 17))
        self.label_2.setObjectName("label_2")
        self.formLayoutWidget = QtWidgets.QWidget(Dialog)
        self.formLayoutWidget.setGeometry(QtCore.QRect(0, 0, 161, 91))
        self.formLayoutWidget.setObjectName("formLayoutWidget")
        self.menuLayout = QtWidgets.QFormLayout(self.formLayoutWidget)
        self.menuLayout.setContentsMargins(0, 0, 0, 0)
        self.menuLayout.setObjectName("menuLayout")
        self.lbinfo = QtWidgets.QLabel(Dialog)
        self.lbinfo.setGeometry(QtCore.QRect(60, 140, 511, 17))
        self.lbinfo.setObjectName("lbinfo")
        self.checkInvert = QtWidgets.QCheckBox(Dialog)
        self.checkInvert.setGeometry(QtCore.QRect(270, 110, 191, 23))
        self.checkInvert.setObjectName("checkInvert")

        self.retranslateUi(Dialog)
        self.listdir.setCurrentRow(2)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "non uniform Magnetic field"))
        self.label.setText(_translate("Dialog", "max B"))
        self.txmaxB.setText(_translate("Dialog", "1000.0"))
        self.btMakeB.setText(_translate("Dialog", "Make B"))
        __sortingEnabled = self.listdir.isSortingEnabled()
        self.listdir.setSortingEnabled(False)
        item = self.listdir.item(0)
        item.setText(_translate("Dialog", "dir X"))
        item = self.listdir.item(1)
        item.setText(_translate("Dialog", "dir Y"))
        item = self.listdir.item(2)
        item.setText(_translate("Dialog", "dir Z"))
        self.listdir.setSortingEnabled(__sortingEnabled)
        self.label_2.setText(_translate("Dialog", "magnetic flow in"))
        self.lbinfo.setText(_translate("Dialog", "info"))
        self.checkInvert.setText(_translate("Dialog", "invert poles (flow)"))

