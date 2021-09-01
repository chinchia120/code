#include <iostream>
using namespace std;
int main(void){
    long long int n;
    while(scanf("%ld",&n)!=EOF){
        if(n<0){
            return 0;
        }else{
            cout<<1+(n*(n+1)/2)<<endl;
        }   
    }
}