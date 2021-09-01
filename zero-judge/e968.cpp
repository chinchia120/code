#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        int arr[3];
        for(int i=0;i<3;i++){
            cin>>arr[i];
        }
        for(int i=n;i>=1;i--){
            if(i==arr[0] || i==arr[1] || i==arr[2]){
                continue;
            }else{
                cout<<"No. "<<i<<endl;
            }
        }
        cout<<"\n";
    }
    return 0;
}