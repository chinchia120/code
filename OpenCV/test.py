import cv2
import numpy as np

# Img_2013
img_2013 = cv2.imread("C:\\Users\\user\\Desktop\\RS\\report2\\F64091091\\2013_supervised_classification_v3.tif")
img_2013 = cv2.resize(img_2013, (300, 300))
img_2013 = img_2013[:, :150]
gray_2013 = cv2.cvtColor(img_2013, cv2.COLOR_BGR2GRAY)
canny_2013 = cv2.Canny(gray_2013, 200, 250)
contours_2013, hierarchy_2013 = cv2.findContours(canny_2013, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

sum_2013 = 0
for cnt in contours_2013:
    area_2013 = cv2.contourArea(cnt)
    sum_2013 += area_2013

print("sum_2013 = {}".format(sum_2013))

cv2.imshow("img_2013", img_2013)
cv2.imshow("gray_2013", gray_2013)
cv2.imshow("canny_2013", canny_2013)

# Img_2016
img_2016 = cv2.imread("C:\\Users\\user\\Desktop\\RS\\report2\\F64091091\\2016_supervised_classification_v3.tif")
img_2016 = cv2.resize(img_2016, (300, 300))
img_2016 = img_2016[:, :150]
gray_2016 = cv2.cvtColor(img_2016, cv2.COLOR_BGR2GRAY)
canny_2016 = cv2.Canny(gray_2016, 200, 250)
contours_2016, hierarchy_2016 = cv2.findContours(canny_2016, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

sum_2016 = 0
for cnt in contours_2016:
    area_2016 = cv2.contourArea(cnt)
    sum_2016 += area_2016

print("sum_2016 = {}".format(sum_2016))

cv2.imshow("img_2016", img_2016)
cv2.imshow("gray_2016", gray_2016)
cv2.imshow("canny_2016", canny_2016)

# Img_2018
img_2018 = cv2.imread("C:\\Users\\user\\Desktop\\RS\\report2\\F64091091\\2018_supervised_classification_v1.tif")
img_2018 = cv2.resize(img_2018, (300, 300))
img_2018 = img_2018[:, :150]
gray_2018 = cv2.cvtColor(img_2018, cv2.COLOR_BGR2GRAY)
canny_2018 = cv2.Canny(gray_2018, 200, 250)
contours_2018, hierarchy_2018 = cv2.findContours(canny_2018, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

sum_2018 = 0
for cnt in contours_2018:
    area_2018 = cv2.contourArea(cnt)
    sum_2018 += area_2018

print("sum_2018 = {}".format(sum_2018))

cv2.imshow("img_2018", img_2018)
cv2.imshow("gray_2018", gray_2018)
cv2.imshow("canny_2018", canny_2018)

# Img_2020
img_2020 = cv2.imread("C:\\Users\\user\\Desktop\\RS\\report2\\F64091091\\2020_supervised_classification_v1.tif")
img_2020 = cv2.resize(img_2020, (300, 300))
img_2020 = img_2020[:, :150]
gray_2020 = cv2.cvtColor(img_2020, cv2.COLOR_BGR2GRAY)
canny_2020 = cv2.Canny(gray_2020, 200, 250)
contours_2020, hierarchy_2020 = cv2.findContours(canny_2020, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

sum_2020 = 0
for cnt in contours_2020:
    area_2020 = cv2.contourArea(cnt)
    sum_2020 += area_2020

print("sum_2020 = {}".format(sum_2020))

cv2.imshow("img_2020", img_2020)
cv2.imshow("gray_2020", gray_2020)
cv2.imshow("canny_2020", canny_2020)

# Show_Img
while True:
        if cv2.waitKey(10) == ord('q'):
            break
