#include <iostream>
using namespace std;
int main(void){
    int n,k=1;
    while(scanf("%d",&n)!=EOF){
        if(n<0){
            return 0;
        }else{
            int tmp=1,cnt=0;
            while(tmp<n){
                tmp*=2;
                cnt++;
            }
            printf("Case %d: %d\n",k,cnt);
            k++;
        }
    }
}