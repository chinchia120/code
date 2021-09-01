#include <iostream>
#include <cmath>
using namespace std;
int main(void){
    double n,m;
    while(scanf("%lf %lf",&n,&m)!=EOF){
        printf("%d\n",(int)log(m)/log(n));
    }
}