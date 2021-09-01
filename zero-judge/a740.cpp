#include <iostream>
using namespace std;
int main(void){
    unsigned long long int ans=0;
    while(cin>>ans){
        int sum=0;
        if(ans==1){
            cout<<"0"<<endl;
            continue;
        }
        for(int i=2;i<ans+1;i++){
            while(ans%i==0){
                ans=ans/i;
                sum=sum+i;
            }
            if(ans==1){
                cout<<sum<<endl;
                break;
            }
        }
    }
    return 0;
}