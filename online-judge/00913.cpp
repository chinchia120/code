#include <iostream>
using namespace std;
int main(void){
    long long int n;
    while(scanf("%ld",&n)!=EOF){
    long long int ans=(3*n*n+6*n+3)/2-9;
    printf("%ld\n",ans);    
    }
    return 0;
}