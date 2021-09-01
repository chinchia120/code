#include <iostream>
using namespace std;

int main(){
    int num;
    while(cin>>num){
        while(num!=0){
            cout<<num<<endl;
            num/=10;
        }
    }
    return 0;
}