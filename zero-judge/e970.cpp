#include <iostream>
using namespace std;
int main(void){
    long long int n;
    while(cin>>n){
        long long int arr[n];
        for(int i=0;i<n;i++){
            cin>>arr[i];
        }
        long long int sum=0;
        for(int i=0;i<n;i++){
            if((i+1)%arr[n-1]==1){
                sum=sum+arr[i];
            }
        }
        //cout<<sum<<endl;
        int q=sum%n;
        if(q==0){
            cout<<n<<" "<<arr[n-1]<<endl;
        }else{
            cout<<q<<" "<<arr[q-1]<<endl;
        }
    }
    return 0;
}