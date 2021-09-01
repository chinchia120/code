#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        int arr[n];
        for(int i=0;i<n;i++){
            cin>>arr[i];
        }
        for(int i=0;i<n;i++){
            for(int j=i+1;j<n;j++){
                    if(arr[i]>arr[j]){
                        int k=0;
                        k=arr[i];
                        arr[i]=arr[j];
                        arr[j]=k;
                    }
            }
        }
        for(int i=0;i<n;i++){
            cout<<arr[i]<<" ";
        }
        cout<<endl;
    }
    return 0;
}