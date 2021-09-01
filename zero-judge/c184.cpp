#include <iostream>
using namespace std;

int main(){
    int ans=0;
    while(cin>>ans){
        if(ans==0){
            return 0;
        }else{
            int sum=0;
            for(int i=1;i<ans;i++){
                if(ans%i==0){
                    sum=sum+i;
                }
            }
            //cout<<sum<<endl;
            if(sum==ans){
                cout<<"="<<ans<<endl;
            }else{
                int total=0;
                for(int i=1;i<sum;i++){
                    if(sum%i==0){
                        total=total+i;
                    }
                }
                //cout<<sum<<endl;
                if(total==ans){
                    cout<<sum<<endl;
                }else{
                    cout<<"0"<<endl;
                }
            }
        }

    }
    return 0;
}