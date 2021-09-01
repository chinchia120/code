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
            int a,b,c;
            cin>>a>>b>>c;
            for(int i=a;i<=a*b*c;i++){
                if(i%a==0 && i%b==0 && i%c==0){
                    cout<<i<<endl;
                    break;
                }
            }
        }
    }
    return 0;
}