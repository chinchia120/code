#include <iostream>
#include <cstdlib>
#include <fstream>
using namespace std;
int main(void){
	// Open Image File
	ifstream ifile("original.bmp",ios::binary);
	
	// Read Header
	unsigned char header[54];
	ifile.read((char*)&header,sizeof(header));		
 	int width=header[18] | header[19]<<8 | header[20]<<16 | header[21]<<24;    
 	int height=header[22] | header[23]<<8 | header[24]<<16 | header[25]<<24; 	
	int bits_px=header[28] | header[29]<<8;
	
	printf("%d*%d*%d\n",width,height,bits_px);	
	
	// Read RGB Data
	unsigned char Original_Data[height][width][bits_px/8];
	for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
            for(int k=0;k<bits_px;k++){
                ifile.read((char*)&Original_Data,sizeof(Original_Data));
            }   
        }         
    }
        
	ifile.close();
	
	// Image Processing
	unsigned char Transform_Data[height][width][bits_px/8];
	for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
            for(int k=0;k<bits_px/8;k++){
                Transform_Data[i][j][k]=Original_Data[i][j][k];
            }	
        }
    }
		
 	int New_width=2*width;
 	header[18]=(New_width>>0) & 0x000000ff;
 	header[19]=(New_width>>8) & 0x000000ff;
 	header[20]=(New_width>>16) & 0x000000ff;
 	header[21]=(New_width>>24) & 0x000000ff;
    
 	int New_height=2*height;
 	header[22]=(New_height>>0) & 0x000000ff;
 	header[23]=(New_height>>8) & 0x000000ff;
 	header[24]=(New_height>>16) & 0x000000ff;
 	header[25]=(New_height>>24) & 0x000000ff;
	 	
	int file_size=New_width*New_height*bits_px/8+54;
 	header[2]=(file_size>>0) & 0x000000ff;
 	header[3]=(file_size>>8) & 0x000000ff;
	header[4]=(file_size>>16) & 0x000000ff;
	header[5]=(file_size>>24) & 0x000000ff; 
	
	printf("%d*%d*%d\n",New_width,New_height,bits_px);	
	
	// Write New Image
	ofstream ofile("Q3.bmp",ios::binary);	
	ofile.write((char*)&header,sizeof(header));
	for(int i=0;i<New_height;i++){
        for(int j=0;j<New_width;j++){
            for(int k=0;k<bits_px/8;k++){
                ofile.write((char*)&Transform_Data[i][j][k],sizeof(Transform_Data[i][j][k]));
                ofile.write((char*)&Transform_Data[i][j][k],sizeof(Transform_Data[i][j][k]));
            }				
        }		
    }
		
	ofile.close();
	
	return 0;
}