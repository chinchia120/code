# -*- coding: utf-8 -*-
import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton)

class MyWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('my window')
        self.setGeometry(500, 500, 400, 400)

        self.mybutton = QPushButton('button', self)
        self.mybutton.move(60, 50)
        self.mybutton.clicked.connect(self.onButtonClick)

    def onButtonClick(self):
        self.mybutton.setText('hello wolrd')

if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MyWidget()
    w.show()
    sys.exit(app.exec_())