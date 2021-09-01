#include <iostream>
using namespace std;

int main(){
    int a,b,d;
    while(cin>>a>>b>>d){
        while(a!=b){
            cout<<a<<" ";
            a=a+d;
        }
        while(a==b){
            cout<<b;
            break;
        }
    }
    return 0;
}