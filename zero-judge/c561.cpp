#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        int total=0;
        for(int i=0;i<n;i++){
            string str;
            cin>>str;
            int dig=str.length(),a=stoi(str);
            int arr[dig];
            for(int j=0;j<dig;j++,a/=10){
                arr[j]=a%10;
                //cout<<arr[j]<<endl;
            }
            int sum=0,k=1;
            for(int j=dig-1;j>=0;j--,k*=10){
                sum=sum+arr[j]*k;
            }
            //cout<<sum<<endl;
            if(sum>total){
                total=sum;
            }
        }
        cout<<total<<endl;
    }
    return 0;
}