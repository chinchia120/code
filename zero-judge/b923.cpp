#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        int stack,num[n],cnt=0;
        for(int i=0;i<n;i++){
            cin>>stack;
            if(stack==1){
                cnt--;
            }else if(stack==2){
                //cout<<cnt<<endl;
                cout<<num[cnt-1]<<endl;
            }else if(stack==3){
                cin>>num[cnt];
                //cout<<num[cnt]<<endl;;
                cnt++;
                //cout<<cnt<<endl;
            }
        }
    }
    return 0;
}