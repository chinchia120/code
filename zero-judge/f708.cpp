#include <iostream>
using namespace std;
int main(void){
    int n,m;
    while(cin>>n>>m){
        int sum1=0,sum2=0,k;
        for(int i=0;i<n;i++){
            cin>>k;
            sum1=sum1+k;
        }
        for(int i=0;i<m;i++){
            cin>>k;
            sum2=sum2+k;
        }
        if(n>m && sum1>sum2){
            cout<<"Yes"<<endl;
        }else{
            cout<<"No"<<endl;
        }
    }
    return 0;
}