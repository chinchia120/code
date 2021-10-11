#include <iostream>
using namespace std;
int main(void){
    int a;
    while(scanf("%d",&a)!=EOF){
        if(a==-99999){
            return 0;
        }
        while(a%17!=0){
            a++;
            if(a%5==0){
                a++;
            }
        }
        printf("%d\n",a);
    }
}