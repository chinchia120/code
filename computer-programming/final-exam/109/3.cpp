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
            int k;
            cin>>k;
            int num[k+2]={};
            for(int j=0;j<k;j++){
                cin>>num[j];
            }
            num[k]=num[k-1];
            num[k+1]=num[k-2];
            int dec=0,inc=0;
            for(int j=0;j<k;j++){
                if(num[j]<num[j+1]){
                    inc++;
                    if(num[j+1]<num[j+2]){
                        inc--;
                    }
                }
                if(num[j]>num[j+1]){
                    dec++;
                    if(num[j+1]>num[j+2]){
                        dec--;
                    }
                }
            }
            cout<<inc<<" "<<dec<<endl;
        }
    }
    return 0;
}
