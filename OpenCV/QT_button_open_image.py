# -*- coding: utf-8 -*-
import sys
import cv2
from PyQt5.QtGui import QImage, QPixmap
from PyQt5.QtWidgets import QApplication, QDialog, QFileDialog, QGridLayout, QLabel, QPushButton

class MyDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.resize(400, 300)
        self.label = QLabel()
        self.btnOpen = QPushButton('Open Image', self)

        layout = QGridLayout(self)
        layout.addWidget(self.label, 0, 0, 4, 4)
        layout.addWidget(self.btnOpen, 4, 0, 1, 1)

        self.btnOpen.clicked.connect(self.openSlot)

    def openSlot(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Open Image', 'Image', '*.png *.jpg *.bmp')
        if filename is '':
            return
        self.img = cv2.imread(filename, -1)
        if self.img.size == 1:
            return
        self.showImage()

    def showImage(self):
        height, width, channel = self.img.shape
        bytesPerline = 3 * width
        self.qImg = QImage(self.img.data, width, height, bytesPerline, QImage.Format_RGB888).rgbSwapped()
        self.label.setPixmap(QPixmap.fromImage(self.qImg))

if __name__ == '__main__':
    a = QApplication(sys.argv)
    dialog = MyDialog()
    dialog.show()
    sys.exit(a.exec_())