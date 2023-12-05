import cv2
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QApplication, QDialog, QFileDialog, QGridLayout, QLabel, QPushButton, QLineEdit

class MyDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Calculate")
        self.resize(1000, 500)
        layout = QGridLayout(self)

        self.label_discribe1 = QLabel('Please select the image', self)
        self.label_discribe1.setGeometry(15, 190, 1000, 60)
        
        self.btn_Load_Image = QPushButton('Load Image', self)
        layout.addWidget(self.btn_Load_Image, 0, 0)
        self.btn_Load_Image.clicked.connect(self.Load_Image)
        
        self.label_Image_Name = QLabel('No image loaded', self)
        self.label_Image_Name.setGeometry(15, 250, 1000, 60)
        
        self.btn_Run = QPushButton('Run', self)
        layout.addWidget(self.btn_Run, 0, 1)
        self.btn_Run.clicked.connect(self.Run)

        self.path = os.path.join(sys.path[0])
        self.imagename = ' '
           
    def Load_Image(self):
        self.imagename, _ = QFileDialog.getOpenFileName(self, "Open Image", self.path)
        #self.label_Image_Name.setText(self.imagename.split('/')[-1])
        self.label_Image_Name.setText(self.imagename)
    
    def Show_Img(self):
        if cv2.waitKey(0) == ord('q'):
            cv2.destroyAllWindows()
      
    def Run(self):
        img = cv2.imread(self.imagename)
        pixel_size = 100*100

        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        thresh, binary = cv2.threshold(gray, 200, 255, cv2.THRESH_BINARY)
        guassian = cv2.GaussianBlur(binary, (7, 7), 0)
        canny = cv2.Canny(guassian, 100, 250)
        contours, hierarchy = cv2.findContours(canny, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

        sum = 0
        for i in contours:
            area = cv2.contourArea(i)
            sum += area

        print("sum  = {:.2f} pixel.".format(sum))
        print("area = {:.2f} km^2.".format(sum*pixel_size))

        cv2.imshow("img", img)
        cv2.imshow("canny", canny)
        self.Show_Img()

if __name__ == '__main__':
    a = QApplication(sys.argv)
    dialog = MyDialog()
    dialog.show()
    sys.exit(a.exec_())