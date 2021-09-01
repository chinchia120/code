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
            string str;
            cin>>str;
            int hr=(str[0]-'0')*10+(str[1]-'0');
            int min=(str[3]-'0')*10+(str[4]-'0');
            //cout<<hr<<" "<<min<<endl;
            float deg=fabs((hr*30+min*0.5)-(min*6));
            if(deg>180){
                printf("%.1f\n",360-deg);
            }else{
                printf("%.1f\n",deg);
            }
        }
    }
    return 0;
}