import cv2

cap = cv2.VideoCapture("C:\\Users\\user\\Documents\\code\\openCV\\video_test.mp4")
# cap = cv2.VideoCapture(0) # 0: notebook camara, 1: external camara

while True:
    ret, frame = cap.read()
    frame = cv2.resize(frame, (0, 0), fx=0.5, fy=0.5)
    if ret:
        cv2.imshow('video', frame)
    else:
        break

    if cv2.waitKey(10) == ord('q'):
        break
