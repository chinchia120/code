#include <iostream>
using namespace std;
int digit(int n){
    int cnt=0;
    while(n!=0){
        n/=10;
        cnt++;
    }
    return cnt;
}
int main(void){
    int n;
    while(cin>>n){
        if(n==0){
            return 0;
        }
        cout<<digit(n)<<endl;
    }
}