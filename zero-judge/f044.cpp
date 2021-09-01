#include <iostream>
using namespace std;
long long int day(long long int k){
    long long int cnt=0;
    while(k!=1){
        cnt++;
        k/=2;
    }
    while(k==1){
        //cout<<cnt<<endl;
        return cnt;
    }
}
int main(void){
    long long int a,b;
    while(cin>>a>>b){
        long long int sum=0;
        if(a==1){
            sum=a+b;
            cout<<day(sum)<<endl;
        }else{
            b=b/a;
            a=1;
            sum=a+b;
            cout<<day(sum)<<endl;
        }
    }
    return 0;
}