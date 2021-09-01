#include <iostream>
#include <algorithm>
using namespace std;
int main(void){
    int n,m;
    while(scanf("%d %d",&n,&m)!=EOF){
        cout<<__gcd(n,m)<<endl;
    }
}