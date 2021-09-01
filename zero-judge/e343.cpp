#include <iostream>
using namespace std;
int main(void){
    float w,h;
    while(cin>>w>>h){
        float BMI=w/h/h*10000;
        //cout<<BMI<<endl;
        printf("%-10.1f\n",BMI);
    }
    return 0;
}