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
    cin>>n;
    int arr[10]={},cnt[10]={};
    for(int i=0;i<n;i++){
        string area;
        int num;
        cin>>area>>num;
        if(area[0]=='A'){
            arr[0]=arr[0]+num;
            cnt[0]++;
        }else if(area[0]=='B'){
            arr[1]=arr[1]+num;
            cnt[1]++;
        }else if(area[0]=='C'){
            arr[2]=arr[2]+num;
            cnt[2]++;
        }else if(area[0]=='D'){
            arr[3]=arr[3]+num;
            cnt[3]++;
        }else if(area[0]=='E'){
            arr[4]=arr[4]+num;
            cnt[4]++;
        }else if(area[0]=='F'){
            arr[5]=arr[5]+num;
            cnt[5]++;
        }else if(area[0]=='G'){
            arr[6]=arr[6]+num;
            cnt[6]++;
        }else if(area[0]=='H'){
            arr[7]=arr[7]+num;
            cnt[7]++;
        }else if(area[0]=='I'){
            arr[8]=arr[8]+num;
            cnt[8]++;
        }else if(area[0]=='J'){
            arr[9]=arr[9]+num;
            cnt[9]++;
        }
    }
    for(int i=0;i<10;i++){
        cout<<(char)(i+65)<<" ";
        if(arr[i]==0){
            cout<<"0 0\n";
        }else{
            cout<<arr[i]<<" ";
            double ave=(float)arr[i]/cnt[i];
            printf("%.1f\n",ave);
        }
    }
    return 0;
}
