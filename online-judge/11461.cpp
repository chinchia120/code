#include <iostream>
#include <cmath>
using namespace std;
int main(void){
    int a,b;
    while(cin>>a>>b){
        if(a==0 && b==0){
            return 0;
        }
        int cnt=0;
        for(int i=a;i<b+1;i++){
            int q=sqrt(i);
            if(i==(q*q)){
                cnt++;
            }
        }
        printf("%d\n",cnt);
    }
}