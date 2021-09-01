#include <iostream>
#include <algorithm>
using namespace std;

int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            int m;
            cin>>m;
            int arr[m];
            cin>>arr[0]>>arr[1]>>arr[2]>>arr[3];
            for(int j=4;j<m;j++){
                arr[j]=arr[j-4]+arr[j-1];
            }
            sort(arr,arr+m);
            cout<<arr[m/2]<<endl;
        }
    }    
}