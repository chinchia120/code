import cv2

def Show_Img(self):
    if cv2.waitKey(0) == ord('q'):
        cv2.destroyAllWindows()

def Count_Contours(img):
    #img = cv2.resize(img, (300, 300))
    img = img[:, :int(img.shape[1]/2)]
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    canny = cv2.Canny(gray, 200, 250)
    contours, hierarchy = cv2.findContours(canny, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    
    sum = 0
    for cnt in contours:
        area = cv2.contourArea(cnt)
        sum += area
    
    print("sum  = {:.2f} pixel.".format(sum))
    print("area = {:.2f} km^2".format(sum * 0.03 * 0.03))
    
    #cv2.imshow("img", img)
    cv2.imshow("canny", canny)
    Show_Img(None)
    
def Resize_Image(path):
    img = cv2.imread(path) 
    img = cv2.resize(img, (300, 300))
    img = img[:, :int(img.shape[1] / 2)]
    
    cv2.imshow('img', img)
    Show_Img(None)

img_2013 = cv2.imread("C:\\Users\\user\\Documents\\code_all\\Remote-Sensing\\report2\\picture\\2013_supervised_classification_v3.tif")
Count_Contours(img_2013)

img_2016 = cv2.imread("C:\\Users\\user\\Documents\\code_all\\Remote-Sensing\\report2\\picture\\2016_supervised_classification_v3.tif")
Count_Contours(img_2016)

img_2018 = cv2.imread("C:\\Users\\user\\Documents\\code_all\\Remote-Sensing\\report2\\picture\\2018_supervised_classification_v1.tif")
Count_Contours(img_2018)

img_2020 = cv2.imread("C:\\Users\\user\\Documents\\code_all\\Remote-Sensing\\report2\\picture\\2020_supervised_classification_v1.tif")
Count_Contours(img_2020)

#path = "C:\\Users\\user\\Documents\\NCKU-Data\\111-1\\Remote-Sensing\\report\\report2\\picture\\ENVI\\2013_NIR_v5.tif"
#Resize_Image(path)

