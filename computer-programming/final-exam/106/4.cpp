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
    while(scanf("%d",&n)!=EOF){
        for(int i=0;i<n;i++){
            int arr[3],a,b;
            for(int j=0;j<3;j++){
                cin>>arr[j];
            }
            cin>>a>>b;
            for(int j=0;j<3;j++){
                for(int k=j+1;k<3;k++){
                    if(arr[j]<arr[k]){
                        int tmp=arr[j];
                        arr[j]=arr[k];
                        arr[k]=tmp;
                    }
                }
            } 
            float score=0;
            if(a>b){
                score=(arr[0]+arr[1])/2*0.6+(a*0.6+b*0.4)*0.4;
                printf("%.2f\n",score);
            }else{
                score=(arr[0]+arr[1])/2*0.6+(b*0.6+a*0.4)*0.4;
                printf("%.2f\n",score);
            }           
        }
    }
    return 0;
}
