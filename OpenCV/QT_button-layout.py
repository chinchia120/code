# -*- coding: utf-8 -*-
import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QGridLayout,
                             QLabel, QLineEdit, QPushButton)

class MyWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('my window')
        self.setGeometry(50, 50, 200, 150)

        layout = QGridLayout()
        self.setLayout(layout)

        self.mylabel = QLabel('Input:', self)
        layout.addWidget(self.mylabel, 0, 0)
        self.mylineedit = QLineEdit(self)
        layout.addWidget(self.mylineedit, 0, 1, 1, 2)
        self.mybutton1 = QPushButton('button 1', self)
        layout.addWidget(self.mybutton1, 1, 0, 1, 3)
        self.mybutton2 = QPushButton('button 2', self)
        layout.addWidget(self.mybutton2, 2, 0, 2, 1)
        self.mybutton3 = QPushButton('button 3', self)
        layout.addWidget(self.mybutton3, 2, 1)
        self.mybutton4 = QPushButton('button 4', self)
        layout.addWidget(self.mybutton4, 2, 2)
        self.mybutton5 = QPushButton('button 5', self)
        layout.addWidget(self.mybutton5, 3, 1)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MyWidget()
    w.show()
    sys.exit(app.exec_())