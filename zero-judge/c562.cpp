#include <iostream>
using namespace std;
int length(long long int k){
    int cnt=0;
    while(k!=0){
        cnt++;
        k=k/10;
    }
    while(k==0){
        return cnt;
    }
}
int main(){
    long long int ans;
    while(cin>>ans){
        if(ans==0){
            cout<<"1"<<endl;
            continue;
        }
        //cout<<ans<<endl;
        int dig=length(ans),cnt=0;
        //cout<<ans<<" "<<dig<<endl;
        int arr[dig];
        for(int i=0;i<dig;i++,ans/=10){
            arr[i]=ans%10;
        }
        for(int i=0;i<dig;i++){
            if(arr[i]==0){
                cnt++;
            }else if(arr[i]==6){
                cnt++;
            }else if(arr[i]==8){
                cnt=cnt+2;
            }else if(arr[i]==9){
                cnt++;
            }
        }
        cout<<cnt<<endl;
    }
    return 0;
}