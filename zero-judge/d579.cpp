#include <iostream>
#include <string>
using namespace std;
int main(void){
    double num;
    while(cin>>num){
        if(num<0){
            printf("|%.4lf|=%.4lf\n",num,-num);
        }else{
            printf("|%.4lf|=%.4lf\n",num,num);
        }
    }
}