#include <iostream>
using namespace std;
int main(void){
    long long int id;
    while(cin>>id){
        int flag=0;
        for(int i=2;i<id;i++){
            if(id%i==0){
                flag=1;
                break;
            }
        }
        if(flag==0){
            cout<<"yes"<<endl;
        }else{
            cout<<"no"<<endl;
        }
    }
    return 0;
}