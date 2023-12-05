import cv2

img = cv2.imread("C:\\Users\\user\\Documents\\code\\openCV\\OpenCV_Logo_with_text.png")
img = cv2.resize(img, (300, 400))
img = cv2.resize(img, (0, 0), fx=0.5, fy=0.5)


