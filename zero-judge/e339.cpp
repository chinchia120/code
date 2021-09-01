#include <iostream>
using namespace std;

int main(void){
    long long int n;
    while(cin>>n){
        long long int sum=0;
        for(long long int i=0;i<n;i++){
            long long int a;
            cin>>a;
            sum=sum+a;
            cout<<sum<<" ";
        }
        cout<<endl;
    }
    return 0;
}