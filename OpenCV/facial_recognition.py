import cv2

# Load Img
img = cv2.imread("C:\\Users\\user\\Documents\\code\\openCV\\lenna.jpg")
#img = cv2.imread("C:\\Users\\user\\Documents\\code\\openCV\\face2.jpg")

# BGR to Gray
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

# Load Facial Recognition Modle
faceCascade = cv2.CascadeClassifier("C:\\Users\\user\\Documents\\code\\openCV\\face_detect.xml")
faceRect = faceCascade.detectMultiScale(gray, 1.1, 7)
print(len(faceRect))

# Creat Rectangle
for (x, y, w, h) in faceRect:
    cv2.rectangle(img, (x, y), (x+w, y+h), (0, 255, 0), 2)

# Img Show
cv2.imshow("img", img)
while True:
    if cv2.waitKey(10) == ord('q'):
        break
