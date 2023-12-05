import cv2
import numpy as np
import random

""" 
img1 = cv2.imread("C:\\Users\\user\\Documents\\code\\openCV\\OpenCV_Logo_with_text.png")
print(type(img1))
print(img1.shape)
"""

""" 
img2 = np.empty((300, 300, 3), np.uint8)

for row in range(300):
    for col in range(300):
        img2[row][col] = [random.randint(0,255), random.randint(0,255), random.randint(0,255)] #BGR

cv2.imshow('img2', img2)
cv2.waitKey(0) 
"""

""" 
img1 = cv2.imread("C:\\Users\\user\\Documents\\code\\openCV\\OpenCV_Logo_with_text.png")


for row in range(300):
    for col in range(img1.shape[1]):
        img1[row][col] = [random.randint(0,255), random.randint(0,255), random.randint(0,255)] #BGR

cv2.imshow('img1', img1)
cv2.waitKey(0)  
"""

img1 = cv2.imread("C:\\Users\\user\\Documents\\code\\openCV\\OpenCV_Logo_with_text.png")

img1_new = img1[:150, :200]

cv2.imshow('img1', img1)
cv2.imshow('img1_new', img1_new)
cv2.waitKey(0)
