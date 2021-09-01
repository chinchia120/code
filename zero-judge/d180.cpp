#include <iostream>
using namespace std;
int main(void){
    int n;
    cin>>n;
    int left=0,right=0;
    for(int i=1;i<=n;i++){
        int num;
        cin>>num;
        if(i%2==1){
            left+=num;
        }else{
            right+=num;
        }
    }
    if(left>right){
        cout<<"left"<<endl;
    }else{
        cout<<"right"<<endl;
    }
}