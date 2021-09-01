#include <iostream>
using namespace std;
int main(void){
    int a,b;
    while(cin>>a>>b){
        if(a==b){
            int c=b-3;
            if((a-c)<c){
                cout<<a-c<<"+"<<c<<"="<<a<<endl;
            }else{
                cout<<c<<"+"<<a-c<<"="<<a<<endl;
            }
        }else{
            if((a-b)<b){
                cout<<a-b<<"+"<<b<<"="<<a<<endl;
            }else{
                cout<<b<<"+"<<a-b<<"="<<a<<endl;
            }
        }
    }
    return 0;
}