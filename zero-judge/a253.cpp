#include <iostream>
using namespace std;
int main(void){
    int no,num,arr[101]={};
    for(int i=0;i<2;i++){ 
        while(cin>>no){
            if(no==-1){
                break;
            }
            cin>>num;
            arr[no]+=num;
        }
    }
    for(int i=0;i<101;i++){
        if(arr[i]==0){
            continue;
        }
        cout<<i<<" "<<arr[i]<<endl;
    }
}