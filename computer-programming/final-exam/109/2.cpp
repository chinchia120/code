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
        float arr[n];
        for(int i=0;i<n;i++){
            float a,b;
            cin>>a>>b;
            arr[i]=arr[i]+a;
            for(int j=0;j<i;j++){
                arr[j]=arr[j]+b/i;
            }
        }
        for(int i=0;i<n;i++){
            printf("%.0f\n",arr[i]);
        }
    }
    return 0;
}