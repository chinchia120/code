import os, sys, cv2
import matplotlib.pyplot as plt

year_Landsat_8 = [2013, 2016, 2018, 2020]
area_Landsat_8 = [375.20, 1096.40, 1913.20, 2890.00]

year_Sentinel_1 = [2016, 2017, 2020]
area_Sentinel_1 = [931.98, 1614.19, 2900.95]

plt.xlabel('year')
plt.ylabel('area(km^2)')
plt.title('Area in Different Year of Landsat 8')
plt.plot(year_Landsat_8, area_Landsat_8)
plt.xticks(year_Landsat_8)
for a,b in zip(year_Landsat_8, area_Landsat_8): 
    plt.text(a, b-100, str(b))
plt.savefig('./picture/Landsat_8.png')
plt.show()

plt.xlabel('year')
plt.ylabel('area(km^2)')
plt.title('Area in Different Year of Sentinel-1')
plt.plot(year_Sentinel_1, area_Sentinel_1)
plt.xticks(year_Sentinel_1)
for a,b in zip(year_Sentinel_1, area_Sentinel_1): 
    plt.text(a, b-100, str(b))
plt.savefig('./picture/Sentinel_1.png')
plt.show()

plt.xlabel('year')
plt.ylabel('area(km^2)')
plt.title('Area in Different Year')
plt.plot(year_Landsat_8, area_Landsat_8, label='Landsat 8', linestyle=":")
#for a,b in zip(year_Landsat_8, area_Landsat_8): 
#    plt.text(a, b+200, str(b))
plt.plot(year_Sentinel_1, area_Sentinel_1, label='Sentinel-1', linestyle="-.")
#for a,b in zip(year_Sentinel_1, area_Sentinel_1): 
#    plt.text(a, b-200, str(b))
plt.legend()
plt.savefig('./picture/plot.png')
plt.show()

def Show_Img(self):
    if cv2.waitKey(0) == ord('q'):
        cv2.destroyAllWindows()

def Count_Contours(img):
    img = cv2.resize(img, (0, 0), fx=0.1, fy=0.1)
        
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    thresh, binary = cv2.threshold(gray, 127, 255, cv2.THRESH_BINARY)
    guassian = cv2.GaussianBlur(binary, (7, 7), 0)
    canny = cv2.Canny(guassian, 200, 250)
    contours, hierarchy = cv2.findContours(canny, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    
    sum = 0
    for i in contours:
        area = cv2.contourArea(i)
        sum += area
    
    print("sum  = {:.2f} pixel.".format(sum))
    print("area = {:.2f} km^2.".format(sum * 0.8 * 0.1))
    print()
    
    cv2.imshow("img", img)
    cv2.imshow("canny", canny)
    Show_Img(None)

path = os.path.join(sys.path[0])

img_2013_ROI = cv2.imread(path + './picture/2013_ROI_v1.png')
Count_Contours(img_2013_ROI)

img_2016_ROI = cv2.imread(path + './picture/2016_ROI_v1.png')
Count_Contours(img_2016_ROI)

img_2018_ROI = cv2.imread(path + './picture/2018_ROI_v1.png')
Count_Contours(img_2018_ROI)

img_2020_ROI = cv2.imread(path + './picture/2020_ROI_v1.png')
Count_Contours(img_2020_ROI)