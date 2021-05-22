# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainWindow123.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog
from PyQt5.QtCore import  QSize
from PyQt5.QtCore import QDir
from PyQt5.QtWidgets import QMessageBox
from SearchnAnalyze import *
from Windows.SearchWindow import Ui_SearchWindow
from Windows.CorrectionWindow import Ui_CorrectionWindow
from Windows.AnalysisWindow import Ui_AnalysisWindow



class Ui_MainWindow(object):

    #Methods for multiple interfaces
    def openSearchWindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_SearchWindow()
        self.ui.setupUi(self.window)
        QtGui.QFontDatabase.addApplicationFont("Resources\GillSansUltraBold.ttf")
        QtGui.QFontDatabase.addApplicationFont("Resources\msuighur.ttf")
        self.window.show()

    def openCorrectionWindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_CorrectionWindow()
        self.ui.setupUi(self.window)
        self.window.show()

    def openAnalysisWindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_AnalysisWindow()
        self.ui.setupUi(self.window)
        self.window.show()
    
    def setupUi(self, MainWindow):
        self.w = 461
        self.h = 496
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(self.w, self.h)
        MainWindow.setFixedSize(self.w, self.h)
        MainWindow.setStyleSheet("background-color:black;")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label1 = QtWidgets.QLabel(self.centralwidget)
        self.label1.setAlignment(QtCore.Qt.AlignCenter)
        self.label1.setGeometry(QtCore.QRect(0, 10, self.w, 20))
        font = QtGui.QFont()
        font.setFamily("Gill Sans Ultra Bold")
        font.setPointSize(12)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label1.setFont(font)
        self.label1.setStyleSheet("color: white;")
        self.label1.setObjectName("label1")
        self.label2 = QtWidgets.QLabel(self.centralwidget)
        self.label2.setGeometry(QtCore.QRect(0, 35, self.w, 20))
        self.label2.setAlignment(QtCore.Qt.AlignCenter)
        font = QtGui.QFont()
        font.setFamily("Gill Sans Ultra Bold")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.label2.setFont(font)
        self.label2.setStyleSheet("color: white;")
        self.label2.setObjectName("label2")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(0, 280, self.w, 20))
        self.label_2.setAlignment(QtCore.Qt.AlignCenter)
        font = QtGui.QFont()
        font.setFamily("Gill Sans Ultra Bold")
        font.setPointSize(11)
        self.label_2.setFont(font)
        self.label_2.setStyleSheet("color: white;")
        self.label_2.setObjectName("label_2")
        self.SearchButton = QtWidgets.QPushButton(self.centralwidget)
        self.SearchButton.setGeometry(QtCore.QRect(self.w//2 - 81//2, 310, 81, 31))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(-1)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.SearchButton.setFont(font)
        self.SearchButton.setStyleSheet("background:blue;\n"
"border-color: rgb(4, 111, 2);\n"
"border-radius:10px;\n"
"font: bold 18px;\n"
"padding: 6px;\n"
"border-style: outset;\n"
"border-width: 0.5px;\n"
"border-color: rgb(0, 0, 127);\n"
"color: white;\n"
"\n"
"")
        self.SearchButton.setObjectName("SearchButton")
        self.CorrectButton = QtWidgets.QPushButton(self.centralwidget)
        self.CorrectButton.setGeometry(QtCore.QRect(self.w//2 - 81//2, 360, 81, 31))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(-1)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.CorrectButton.setFont(font)
        self.CorrectButton.setStyleSheet("background:blue;\n"
"border-color: rgb(4, 111, 2);\n"
"border-radius:10px;\n"
"font: bold 18px;\n"
"padding: 6px;\n"
"border-style: outset;\n"
"border-width: 0.5px;\n"
"border-color: rgb(0, 0, 127);\n"
"color: white;\n"
"")
        self.CorrectButton.setObjectName("CorrectButton")
        self.AnalyseButton = QtWidgets.QPushButton(self.centralwidget)
        self.AnalyseButton.setGeometry(QtCore.QRect(self.w//2 - 81//2, 410, 81, 31))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(-1)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.AnalyseButton.setFont(font)
        self.AnalyseButton.setStyleSheet("background:blue;\n"
"border-color: rgb(4, 111, 2);\n"
"border-radius:10px;\n"
"font: bold 17px;\n"
"padding: 6px;\n"
"border-style: outset;\n"
"border-width: 0.5px;\n"
"border-color: rgb(0, 0, 127);\n"
"color: white;\n"
"")
        self.AnalyseButton.setObjectName("AnalyseButton")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(0, 70, self.w, 181))
        self.label.setText("")
        self.label.setPixmap(QtGui.QPixmap("Resources/dna.jpg"))
        self.label.setScaledContents(True)
        self.label.setObjectName("label")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, self.w, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        #Connect forms
        self.SearchButton.clicked.connect(self.openSearchWindow)
        self.SearchButton.clicked.connect(MainWindow.close)
        self.CorrectButton.clicked.connect(self.openCorrectionWindow)
        self.CorrectButton.clicked.connect(MainWindow.close)
        self.AnalyseButton.clicked.connect(self.openAnalysisWindow)
        self.AnalyseButton.clicked.connect(MainWindow.close)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label1.setText(_translate("MainWindow", "DNA Bases Search Engine,"))
        self.label2.setText(_translate("MainWindow", " Correction and Analysis"))
        self.label_2.setText(_translate("MainWindow", "Click to perform an operation:"))
        self.SearchButton.setText(_translate("MainWindow", "Search"))
        self.CorrectButton.setText(_translate("MainWindow", "Correct"))
        self.AnalyseButton.setText(_translate("MainWindow", "Analyse"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
