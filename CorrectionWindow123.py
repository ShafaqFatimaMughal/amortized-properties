# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'CorrectionWindow123.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog
from PyQt5.QtCore import QDir
from PyQt5.QtWidgets import QMessageBox
from Implementation import *


class Ui_CorrectionWindow(object):

    def openMainWindow(self):
        from MainWindow123 import Ui_MainWindow
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.window)
        self.window.show()

    def setupUi(self, CorrectionWindow):
        CorrectionWindow.setObjectName("CorrectionWindow")
        CorrectionWindow.resize(666, 556)
        CorrectionWindow.setFixedSize(666, 556)
        CorrectionWindow.setStyleSheet("\n"
"background: black;\n"
"")
        self.centralwidget = QtWidgets.QWidget(CorrectionWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label5_2 = QtWidgets.QLabel(self.centralwidget)
        self.label5_2.setGeometry(QtCore.QRect(340, 160, 31, 21))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(14)
        self.label5_2.setFont(font)
        self.label5_2.setStyleSheet("color: white;")
        self.label5_2.setObjectName("label5_2")
        self.label6_2 = QtWidgets.QLabel(self.centralwidget)
        self.label6_2.setGeometry(QtCore.QRect(40, 240, 271, 16))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(14)
        self.label6_2.setFont(font)
        self.label6_2.setStyleSheet("color: white;")
        self.label6_2.setObjectName("label6_2")
        self.textEdit5_2 = QtWidgets.QTextEdit(self.centralwidget)
        self.textEdit5_2.setGeometry(QtCore.QRect(310, 240, 141, 21))
        self.textEdit5_2.setStyleSheet("color: white;")
        self.textEdit5_2.setObjectName("textEdit5_2")
        self.label2_2 = QtWidgets.QLabel(self.centralwidget)
        self.label2_2.setGeometry(QtCore.QRect(40, 110, 461, 21))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(18)
        self.label2_2.setFont(font)
        self.label2_2.setStyleSheet("color: white;")
        self.label2_2.setObjectName("label2_2")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(40, 160, 51, 20))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(14)
        self.label_2.setFont(font)
        self.label_2.setStyleSheet("color: white;")
        self.label_2.setObjectName("label_2")
        self.CorrectButton_2 = QtWidgets.QPushButton(self.centralwidget)
        self.CorrectButton_2.setGeometry(QtCore.QRect(300, 290, 101, 31))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(-1)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.CorrectButton_2.setFont(font)
        self.CorrectButton_2.setStyleSheet("background:blue;\n"
"border-color: rgb(4, 111, 2);\n"
"border-radius:10px;\n"
"font: bold 18px;\n"
"padding: 6px;\n"
"border-style: outset;\n"
"border-width: 0.3px;\n"
"border-color: rgb(0, 0, 127);\n"
"color: white;\n"
"\n"
"")
        self.CorrectButton_2.setObjectName("CorrectButton_2")
        self.textEdit3_2 = QtWidgets.QTextEdit(self.centralwidget)
        self.textEdit3_2.setGeometry(QtCore.QRect(90, 150, 231, 51))
        self.textEdit3_2.setStyleSheet("color: white;")
        self.textEdit3_2.setObjectName("textEdit3_2")
        self.textEdit4_2 = QtWidgets.QTextEdit(self.centralwidget)
        self.textEdit4_2.setGeometry(QtCore.QRect(380, 150, 241, 51))
        self.textEdit4_2.setStyleSheet("color: white;")
        self.textEdit4_2.setObjectName("textEdit4_2")
        self.textEdit_2 = QtWidgets.QTextEdit(self.centralwidget)
        self.textEdit_2.setGeometry(QtCore.QRect(40, 340, 581, 71))
        self.textEdit_2.setStyleSheet("color: white;")
        self.textEdit_2.setObjectName("textEdit_2")
        self.label1_2 = QtWidgets.QLabel(self.centralwidget)
        self.label1_2.setGeometry(QtCore.QRect(220, 10, 241, 31))
        font = QtGui.QFont()
        font.setFamily("Gill Sans Ultra Bold")
        font.setPointSize(18)
        self.label1_2.setFont(font)
        self.label1_2.setStyleSheet("color: white;")
        self.label1_2.setObjectName("label1_2")
        self.GoBack = QtWidgets.QPushButton(self.centralwidget)
        self.GoBack.setGeometry(QtCore.QRect(40, 470, 75, 23))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(-1)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.GoBack.setFont(font)
        self.GoBack.setStyleSheet("background:blue;\n"
"border-color: rgb(4, 111, 2);\n"
"border-radius:10px;\n"
"font: bold 18px;\n"
"padding: 6px;\n"
"border-style: outset;\n"
"border-width: 0.3px;\n"
"border-color: rgb(0, 0, 127);\n"
"color: white;\n"
"\n"
"")
        self.GoBack.setObjectName("GoBack")
        self.DataSetButton = QtWidgets.QPushButton(self.centralwidget)
        self.DataSetButton.setGeometry(QtCore.QRect(40, 80, 101, 21))
        font = QtGui.QFont()
        font.setFamily("Microsoft Uighur")
        font.setPointSize(-1)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.DataSetButton.setFont(font)
        self.DataSetButton.setStyleSheet("background:blue;\n"
"border-color: rgb(4, 111, 2);\n"
"border-radius:10px;\n"
"font: bold 18px;\n"
"padding: 6px;\n"
"border-style: outset;\n"
"border-width: 0.3px;\n"
"border-color: rgb(0, 0, 127);\n"
"color: white;\n"
"\n"
"")
        self.DataSetButton.setObjectName("DataSetButton")
        self.textEdit_3 = QtWidgets.QTextEdit(self.centralwidget)
        self.textEdit_3.setGeometry(QtCore.QRect(150, 80, 241, 21))
        self.textEdit_3.setStyleSheet("color: white;")
        self.textEdit_3.setObjectName("textEdit_3")
        CorrectionWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(CorrectionWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 666, 21))
        self.menubar.setObjectName("menubar")
        CorrectionWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(CorrectionWindow)
        self.statusbar.setObjectName("statusbar")
        CorrectionWindow.setStatusBar(self.statusbar)

        self.retranslateUi(CorrectionWindow)
        QtCore.QMetaObject.connectSlotsByName(CorrectionWindow)

        self.textEdit_3.setDisabled(True)
        #connect button to function
        #self.CorrectButton_2.clicked.connect(lambda: self.Correct(self.textEdit3_2.toPlainText(), self.textEdit4_2.toPlainText(), int(self.textEdit5_2.toPlainText())))
        self.CorrectButton_2.clicked.connect(self.cl)
        self.GoBack.clicked.connect(self.openMainWindow)
        self.GoBack.clicked.connect(CorrectionWindow.close)
        self.DataSetButton.clicked.connect(self.openDialogBox)

    def retranslateUi(self, CorrectionWindow):
        _translate = QtCore.QCoreApplication.translate
        CorrectionWindow.setWindowTitle(_translate("CorrectionWindow", "CorrectionWindow"))
        self.label5_2.setText(_translate("CorrectionWindow", "TO"))
        self.label6_2.setStatusTip(_translate("CorrectionWindow", "As to avoid changing the sequences at many occurences"))
        self.label6_2.setText(_translate("CorrectionWindow", "Enter an Apprioximate Position of this sequence:"))
        self.textEdit5_2.setStatusTip(_translate("CorrectionWindow", "As to avoid changing the sequences at many occurences"))
        self.label2_2.setText(_translate("CorrectionWindow", "Correct bases sequence:"))
        self.label_2.setText(_translate("CorrectionWindow", "FROM"))
        self.CorrectButton_2.setText(_translate("CorrectionWindow", "Correct"))
        self.textEdit3_2.setHtml(_translate("CorrectionWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Sequence1</p></body></html>"))
        self.textEdit4_2.setHtml(_translate("CorrectionWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Sequence2</p></body></html>"))
        self.textEdit_2.setHtml(_translate("CorrectionWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Results</p></body></html>"))
        self.label1_2.setText(_translate("CorrectionWindow", "Gene Correction"))
        self.GoBack.setText(_translate("CorrectionWindow", "Go back"))
        self.DataSetButton.setText(_translate("CorrectionWindow", "Choose a dataset"))


    def openDialogBox(self):
        filename = QFileDialog.getOpenFileName()
        path = filename[0]
        self.textEdit_3.setText(str(path))

    def Correct(self, string1, string2, approxpos=None):
        if self.textEdit_3.toPlainText()=="":
             self.MessageBoxDirectoryError()

        elif self.textEdit3_2.toPlainText()=="" or self.textEdit3_2.toPlainText()=="Sequence1" or self.textEdit4_2.toPlainText()=="" or self.textEdit4_2.toPlainText()=="Sequence2" or ((self.textEdit3_2.toPlainText()=="" and self.textEdit4_2.toPlainText()=="")) or ((self.textEdit3_2.toPlainText()=="Sequence1" and self.textEdit4_2.toPlainText()=="Sequence2")):
             self.MessageBoxGeneError()

        else:
             self.Imp = Implementation(self.textEdit_3.toPlainText())
             self.textEdit_2.setText(str((self.Imp.correction(string1, string2, approxpos))))

    def cl(self):
        try:
            self.Correct(self.textEdit3_2.toPlainText(), self.textEdit4_2.toPlainText(), int(self.textEdit5_2.toPlainText()))
        except:
            self.Correct(self.textEdit3_2.toPlainText(), self.textEdit4_2.toPlainText())
    
    def MessageBoxDirectoryError(self):
        msg = QMessageBox()
        msg.setText("Please choose a data set!")
        msg.setWindowTitle("Error!")
        msg.setIcon(QMessageBox.Critical)
        msg.setStyleSheet("background-color: rgb(0, 0, 0); color: rgb(255, 255, 255)")
        msg.exec_()

    def MessageBoxGeneError(self):
        msg = QMessageBox()
        msg.setText("Cannot proceed until DNA sequence given!")
        msg.setWindowTitle("Error!")
        msg.setIcon(QMessageBox.Critical)
        msg.setStyleSheet("background-color: rgb(0, 0, 0); color: rgb(255, 255, 255)")
        msg.exec_()

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    CorrectionWindow = QtWidgets.QMainWindow()
    ui = Ui_CorrectionWindow()
    ui.setupUi(CorrectionWindow)
    CorrectionWindow.show()
    sys.exit(app.exec_())
