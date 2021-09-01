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
            int a,b,sum=0;
            cin>>a>>b;
            for(int i=a;i<b+1;i++){
                for(int j=i+1;j<b+1;j++){
                    int flag=1;
                    for(int k=2;k<=j;k++){
                        if(i%k==0 && j%k==0){
                            flag=0;
                        }
                    }
                    if(flag==1){
                        sum++;
                        //cout<<i<<" "<<j<<endl;
                    }
                }
            }
            cout<<sum<<endl;
        }
    }
    return 0;
}