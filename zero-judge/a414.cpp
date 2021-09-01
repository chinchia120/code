#include <iostream>
using namespace std;
int main(void){
    long long int ans=0;
    while(scanf("%d",ans)){
        int cnt=0;
        if(ans==0){
            return 0;
        }else{
            while(ans!=1){
                ans=ans/2;
                cnt++;
            }
            if(ans==1){
                printf("%d/n",cnt);
                break;
            }
        }
    }
}