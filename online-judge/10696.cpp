#include <iostream>
using namespace std;
int f91(int n){
    if(n<=100){
        return f91(f91(n+11));
    }else{
        return n-10;
    }
}
int main(void){
    long long int a;
    while(cin>>a){
        if(a==0){
            return 0;
        }else{
            printf("f91(%d) = %d\n",a,f91(a));
        }
    }
}