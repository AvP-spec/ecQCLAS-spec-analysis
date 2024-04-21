# C:\myProgramms\Python\Python312\python test.py
# D:\Andrei_home\python\avp-las01\avp-las01\venv\Scripts\python test.py
# set PATH=%PATH%;C:\your\path\here\
from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QApplication
import sys

app = QApplication(sys.argv)
dname = QtWidgets.QFileDialog.getExistingDirectory()

# dname = QtWidgets.QFileDialog.getOpenFileName(parent=self, options=QtWidgets.QFileDialog.DontUseNativeDialog)
# dname = QtWidgets.QFileDialog.getOpenFileName(options=QtWidgets.QFileDialog.DontUseNativeDialog)
print(dname)