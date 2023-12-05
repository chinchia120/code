import cv2
import numpy as np

img = cv2.imread("C:\\Users\\user\\Documents\\code\\openCV\\cat.jpg")
img = cv2.resize(img, (0, 0), fx=0.5, fy=0.5)

# BGR to Gray
img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

# Gaussian Blur
img_blur = cv2.GaussianBlur(img, (5, 5), 10) # img_name, kernal(odd), standard deviation

# Canny
img_canny = cv2.Canny(img, 100, 200) # img_name, lower_boundary, upper_boundary

# Dilate
kernal_dilate = np.ones((5, 5), np.uint8)
img_dilate = cv2.dilate(img_canny, kernal_dilate, iterations=1)

# Erode
kernal_erode = np.ones((5, 5), np.uint8)
img_erode = cv2.erode(img_dilate, kernal_erode, iterations=1)

# Img Show
cv2.imshow("img", img)
cv2.imshow("img_gray", img_gray)
cv2.imshow("img_blur", img_blur)
cv2.imshow("img_canny", img_canny)
cv2.imshow("img_dilate", img_dilate)
cv2.imshow("img_erode", img_erode)

while True:
    if cv2.waitKey(10) == ord('q'):
        break
