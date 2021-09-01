#include <iostream>
using namespace std;
int main(void){
    int ans=0;
    while(cin>>ans){
        for(int i=2;i<ans+1;i++){
            int cnt=0;
            while(ans%i==0){
                ans/=i;
                cnt++;
            }
            if(cnt==0){
                continue;
            }else if(cnt==1){
                cout<<i;
                if(ans!=1){
                    cout<<" * ";
                }else if(ans==1){
                    break;
                }
            }else{
                cout<<i<<"^"<<cnt;
                if(ans!=1){
                    cout<<" * ";
                }else if(ans==1){
                    break;
                }
            }
        }
    }
    return 0;
}