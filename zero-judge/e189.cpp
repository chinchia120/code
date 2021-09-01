#include <iostream>
using namespace std;
int main(void){
    long long int ans;
    while(cin>>ans){
        if(ans==0){
            cout<<"YES\n";
            continue;
        }
        if(ans%3==0){
            cout<<"YES\n";
        }else{
            cout<<"NO\n";
        }
    }
    return 0;
}