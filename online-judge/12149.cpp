#include <iostream>
using namespace std;
int main(void){
    long long int a;
    while(cin>>a){
        if(a==0){
            return 0;
        }else{
            cout<<(a)*(a+1)*(2*a+1)/6<<endl;
        }
    }
}