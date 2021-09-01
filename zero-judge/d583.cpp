#include <iostream>
using namespace std;

int main(){
    int a,t;
    while(cin>>a){
        for(int i=1;i<a+1;i++){
            cin>>t;
            cout<<i<<" ";
        }
        cout<<endl;
    }
    return 0;
}