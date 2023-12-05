import os
import cv2

def Show_Img(self):
    if cv2.waitKey(0) == ord('q'):
        cv2.destroyAllWindows()

def Count_Contours(img):
    img = cv2.resize(img, (300, 300))
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    canny = cv2.Canny(gray, 200, 250)
    contours, hierarchy = cv2.findContours(canny, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    
    sum = 0
    for cnt in contours:
        area = cv2.contourArea(cnt)
        sum += area
    
    print("sum = {} pixel.".format(sum))
    
    cv2.imshow("img", img)
    cv2.imshow("canny", canny)
    Show_Img(None)

img_2017_ascending = cv2.imread("C:\\Users\\user\\Documents\\code_all\\Remote-Sensing\\report3\\picture\\2016_ascending_v1.tif")
Count_Contours(img_2017_ascending)

img_2017_ascending = cv2.imread("C:\\Users\\user\\Documents\\code_all\\Remote-Sensing\\report3\\picture\\2017_ascending_v1.tif")
Count_Contours(img_2017_ascending)

img_2020_ascending = cv2.imread("C:\\Users\\user\\Documents\\code_all\\Remote-Sensing\\report3\\picture\\2020_ascending_v1.tif")
Count_Contours(img_2020_ascending)
