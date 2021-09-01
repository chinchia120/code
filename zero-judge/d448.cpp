#include <iostream>
using namespace std;
int main(void){
    double t1,t2,t3,x1,x3;
    while(scanf("%lf %lf %lf %lf %lf",&t1,&t2,&t3,&x1,&x3)!=EOF){
        printf("%.6lf\n",(t2-t3)*(x1-x3)/(t1-t3)+x3);
    }
}