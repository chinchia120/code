#include <iostream>
using namespace std;
int main(void){
    while(1){
        float arr[4];
        for(int i=0;i<4;i++){
            cin>>arr[i];
            //cout<<arr[i]<<" ";
        }
        if(arr[0]==0 && arr[1]==0 && arr[2]==0 && arr[3]==0){
            return 0;
        }else{
            float det=(arr[0]*arr[3])-(arr[1]*arr[2]);
            //cout<<det<<endl;
            if(det==0){
                cout<<"cheat!"<<endl;
            }else{
                float rev[4];
                rev[0]=1/det*arr[3];
                rev[1]=1/det*(-arr[1]);
                rev[2]=1/det*(-arr[2]);
                rev[3]=1/det*arr[0];
                printf("%.5f %.5f\n",rev[0],rev[1]);
                printf("%.5f %.5f\n",rev[2],rev[3]);
            }
        }
    }
}
