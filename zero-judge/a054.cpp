#include <iostream>
using namespace std;
int main(void){
    int id;
    while(cin>>id){
        int arr[9];
        for(int i=8;i>=0;i--,id/=10){
            arr[i]=id%10;
        }
        int k=8,sum=0;
        for(int i=0;i<8;i++,k--){
            //cout<<arr[i]<<" ";
            sum=sum+arr[i]*k;
        }
        sum=sum+arr[8];
        //cout<<sum<<endl;
        if((sum+10)%10==0){
            cout<<"BNZ"<<endl;
        }else if((sum+1)%10==0){
            cout<<"AMW"<<endl;
        }else if((sum+2)%10==0){
            cout<<"KLY"<<endl;
        }else if((sum+3)%10==0){
            cout<<"JVX"<<endl;
        }else if((sum+4)%10==0){
            cout<<"HU"<<endl;
        }else if((sum+5)%10==0){
            cout<<"GT"<<endl;
        }else if((sum+6)%10==0){
            cout<<"FS"<<endl;
        }else if((sum+7)%10==0){
            cout<<"ER"<<endl;
        }else if((sum+8)%10==0){
            cout<<"DOQ"<<endl;
        }else if((sum+9)%10==0){
            cout<<"CIP"<<endl;
        }

    }
    return 0;
}