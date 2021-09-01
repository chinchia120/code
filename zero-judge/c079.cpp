#include <iostream>
using namespace std;
int main(void){
    int a,b;
    while(cin>>a>>b){
        cout<<a+a/b+(a/b+a%b)/b<<endl;
    }
    return 0;
}