#include <iostream>
using namespace std;
int main(void){
    int id;
    while(cin>>id){
        int arr[9];
        for(int i=0;i<9;i++,id/=10){
            arr[i]=id%10;
        }
        int max[2]={0,0};
        for(int i=0;i<9;i++){
            if(arr[i]>max[0]){
                max[0]=arr[i];
            }else if(arr[i]>max[1] && arr[i]<=max[0]){
                max[1]=arr[i];
            }
        }
        //cout<<max[0]<<" "<<max[1]<<endl;
        int sum=max[0]*max[0]+max[1]*max[1];
        //cout<<sum<<endl;
        int check=arr[0]+arr[1]*10+arr[2]*100;
        //cout<<check<<endl;
        if(check==sum){
            cout<<"Good Morning!\n";
        }else{
            cout<<"SPY!\n";
        }
    }
    return 0;
}