#include <iostream>
#include <math.h>
using namespace std;
int main(void){
    long long unsigned int n;
    while(scanf("%d",&n)){
        int flag=0;
        for(long long unsigned int i=2;i<n;i++){
            if(n%i==0){
                flag=1;
                break;
            }
        }
        if(flag==1){
            printf("非質數\n");
        }else{
            printf("質數\n");
        }
    }
    return 0;
}