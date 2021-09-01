#include <iostream>
using namespace std;
long long int F(int n){
    if(n==1){
        return 1;
    }else if(n==2){
        return 2;
    }else{
        return F(n-1)+F(n-2);
    }
}
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        if(n==0){
            return 0;
        }
        printf("%ld\n",F(n));
    }
}
