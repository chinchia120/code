import cv2
import numpy as np

# Generate a Black Photo
img = np.zeros((600, 600, 3), np.uint8)

# Line
cv2.line(img, (0, 0), (400, 300), (255, 0, 0), 1)
cv2.line(img, (0, 0), (img.shape[1], img.shape[0]), (0, 255, 0), 2)

# Rectangle
cv2.rectangle(img, (0, 0), (400, 300), (0, 0, 255), cv2.FILLED)

# Circle
cv2.circle(img, (300, 400), 30, (255, 0, 0), 4)

# Text
cv2.putText(img, "Hello", (100, 500), cv2.FONT_HERSHEY_SCRIPT_SIMPLEX, 2, (255, 255, 255), 1) # doesn't support Chinese

# Img Show
cv2.imshow("img", img)
while True:
    if cv2.waitKey(10) == ord('q'):
        break
