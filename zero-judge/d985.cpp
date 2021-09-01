#include <iostream>
using namespace std;
int main(void){
    int track;
    while(cin>>track){
        for(int i=1;i<=track;i++){
            int n;
            cin>>n;
            int arr[n][2],sum=0,small[n]={};
            for(int j=0;j<n;j++){
                for(int k=0;k<2;k++){
                    cin>>arr[j][k];
                    //cout<<arr[j][k]<<endl;
                    if(k==0){
                        sum=sum+arr[j][0]*60;
                        small[j]=small[j]+arr[j][0]*60;
                    }else if(k==1){
                        sum=sum+arr[j][1];
                        small[j]=small[j]+arr[j][1];
                    }
                }
                //cout<<small[j+1]<<endl;
            }
            //cout<<sum<<endl;
            int best=small[0];
            for(int j=0;j<n;j++){
                if(best>small[j]){
                    best=small[j];
                }
            }
            int average=sum/n;
            printf("Track %d:\n",i);
            printf("Best Lap: %d minute(s) %d second(s).\n",best/60,best%60);
            printf("Average: %d minute(s) %d second(s).\n",average/60,average%60);
        }    
    }
    return 0;
}