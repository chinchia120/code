#include <iostream>
using namespace std;
int main(void){
    double H,M,S;
    while(cin>>H>>M>>S){
        if((M/S)<=H){
            cout<<(int)H<<" "<<(int)M<<" "<<(int)S<<". I will make it!\n";
        }else{
            cout<<(int)H<<" "<<(int)M<<" "<<(int)S<<". I will be late!\n";
        }
    }
    return 0;
}