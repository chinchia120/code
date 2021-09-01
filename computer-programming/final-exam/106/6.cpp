#include <bits\stdc++.h>
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

int main() {
    int n;
    while(scanf("%d",n)!=EOF){
        for(int i=0;i<n;i++){
            int arr[5];
            for(int j=0;j<5;j++){
                cin>>arr[j];
            }
            for(int j=0;j<5;j++){
                for(int k=j+1;k<5;k++){
                    if(arr[j]>arr[k]){
                        int tmp=arr[j];
                        arr[j]=arr[k];
                        arr[k]=tmp;
                    }
                }
            }
            if(arr[0]==1 && arr[0]==10 && arr[0]==11 && arr[0]==12 && arr[0]==13){
                cout<<"Yes"<<endl;
            }else{
                int flag=1;
                for(int j=0;j<4;j++){
                    if(arr[j+1]-arr[j]!=1){
                        flag=0;
                        break;
                    }
                }
                if(flag==1){
                    cout<<"Yes"<<endl;
                }else{
                    cout<<"No"<<endl;
                }
            }
        }
    }
    return 0;
}
