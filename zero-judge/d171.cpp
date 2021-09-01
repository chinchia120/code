#include <iostream>
#include <cmath>
using namespace std;
int main(void){
    int a,b;
    while(cin>>a>>b){
        cout<<(int)(b*log10(a))+1<<endl;
    }
}