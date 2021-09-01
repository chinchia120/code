#include <iostream>
using namespace std;

int main(){
    int a,b;
    while(cin>>a>>b){
        if(b>=a){
            cout<<b-a<<endl;
        }else if(a>b){
            cout<<100-(a-b)<<endl;
        }
    }
    return 0;
}