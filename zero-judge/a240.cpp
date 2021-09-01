#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            int dig;
            cin>>dig;
            int sum=0,q=1,tmp;
            for(int i=0;i<dig;i++){
                q*=10;
                tmp=q/17;
                sum=sum+tmp;
                q%=17;
            }
            cout<<tmp<<" "<<sum<<endl;
        }
    }
    return 0;
}
