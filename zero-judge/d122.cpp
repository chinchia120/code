#include <iostream>
using namespace std;
int count(long int n){
    int cnt=0,tmp=n;
    for(int i=5;i<=tmp;i=i*5){
        cnt=cnt+n/i;
    }
    return cnt;
}
int main(void){
    long int n;
    while(scanf("%ld",&n)!=EOF){
        printf("%d\n",count(n));
    }
    return 0;
}