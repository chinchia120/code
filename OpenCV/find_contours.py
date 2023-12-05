import cv2

img = cv2.imread("C:\\Users\\user\\Documents\\code\\openCV\\shape.png")
#img = cv2.imread("C:\\Users\\user\\Documents\\code\\openCV\\shape2.jpg")
img = cv2.resize(img, (0, 0), fx=2, fy=2)

img_Contour = img.copy()

img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
canny = cv2.Canny(img, 150, 200)
contours, hierarchy = cv2.findContours(canny, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

for cnt in contours:
    cv2.drawContours(img_Contour, cnt, -1, (255, 0, 0), 1)
    
    area = cv2.contourArea(cnt)
      
    if area > 0:
        peri = cv2.arcLength(cnt, True)
        vertices = cv2.approxPolyDP(cnt, peri * 0.02, True)
        corners = len(vertices)
        x, y, w, h = cv2.boundingRect(vertices)
        cv2.rectangle(img_Contour, (x, y), (x+w, y+h), (0, 255, 0), 4)
        
        if corners == 3:
            cv2.putText(img_Contour, "triangle", (x, y-5), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)
        elif corners == 4:
            cv2.putText(img_Contour, "rectangle", (x, y-5), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)
        elif corners == 5:
            cv2.putText(img_Contour, "pentagon", (x, y-5), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)
        elif corners == 6:
            cv2.putText(img_Contour, "hexagon", (x, y-5), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)
        elif corners > 6:
            cv2.putText(img_Contour, "circle", (x, y-5), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)


cv2.imshow("img", img)
cv2.imshow("canny", canny)
cv2.imshow("img_Contour", img_Contour)

while True:
    if cv2.waitKey(10) == ord('q'):
        break
