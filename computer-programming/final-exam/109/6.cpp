#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <stack>
#include <algorithm>
#include <vector>
using namespace std;
bool check(int arr[3][3]){
    if((arr[0][0]==0 && arr[0][1]==0 && arr[0][2]==0) ||
    (arr[1][0]==0 && arr[1][1]==0 && arr[1][2]==0) ||
    (arr[2][0]==0 && arr[1][1]==0 && arr[2][2]==0) ||
    (arr[0][0]==0 && arr[1][0]==0 && arr[2][0]==0) ||
    (arr[0][1]==0 && arr[1][1]==0 && arr[2][1]==0) ||
    (arr[0][2]==0 && arr[1][2]==0 && arr[2][2]==0) ||
    (arr[0][0]==0 && arr[1][1]==0 && arr[2][2]==0) ||
    (arr[0][2]==0 && arr[1][1]==0 && arr[2][0]==0)){
    	return true;
	}else{
		return false;
	}
}
int main() {
    int n;
    cin>>n;
    for(int i=0;i<n;i++){
        string str[3];
        int arr[3][3];
        for(int j=0;j<3;j++){
            cin>>str[j];
            for(int k=0;k<3;k++){
                arr[j][k]=str[j].at(k)-'0';
            }
        }
        int flag=0;
        for(int j=1;j<10;j++){
            for(int k=0;k<3;k++){
                for(int l=0;l<3;l++){
                    if(j==arr[k][l]){
                        arr[k][l]=0;
                        if(check(arr)){
                            cout<<j<<endl;
                            flag=1;
                            break;
                        }
                    }
                }
                if(flag==1){
                    break;
                }
            }
            if(flag==1){
                break;
            }
        }
    }
    return 0;
}
