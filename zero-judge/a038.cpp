#include <iostream>
using namespace std;

int main(){
    long long int a;
    while(cin>>a){
        while(a%10==0){
            a=a/10;
            if(a==0){
                cout<<"0"<<endl;
                break;
            }
        }
        while(a>0){
            cout<<a%10;
            a=a/10;
        }
        cout<<endl;
    }
    return 0;
}