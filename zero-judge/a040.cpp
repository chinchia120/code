#include <iostream>
#include <math.h>
using namespace std;
int dig(int k){
    int cnt=0;
    while(k!=0){
        cnt++;
        k/=10;
    }
    while(k==0){
        return cnt;
    }
}
int main(void) {
    int a,b;
    while(cin>>a>>b){
    	int flag=0;
        for(int i=a;i<=b;i++){
            int tmp=i,length=dig(tmp);
            //cout<<length<<endl;
            int arr[length];
            for(int j=0;j<length;j++,tmp/=10){
                arr[j]=tmp%10;
                //cout<<arr[j]<<" ";
            }
            //cout<<"\n";
            int sum[length];
            for(int j=0;j<length;j++){
                sum[j]=pow(arr[j],length);
                //cout<<sum[j]<<" ";
            }
            int total=0;
            for(int j=0;j<length;j++){
                total=total+sum[j]; 
            }
            //cout<<total<<endl;
            if(total==i){
                cout<<i<<" ";
                flag=1;
            }
        }
        if(flag==0){
            cout<<"none";
        }
        cout<<endl;
    }
    return 0;
}