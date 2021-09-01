#include <iostream>
using namespace std;
int dig(int k){
    int cnt=0;
    while(k!=0){
        k=k/10;
        cnt++;
    }
    while(k==0){
        return cnt;
    }
}
int main(void){
    int a;
    while(cin>>a){
        int length=dig(a);
        cout<<length<<endl;
        int arr[length];
        for(int i=0;i<length;i++,a/=10){
            arr[i]=a%10;
            cout<<arr[i]<<" ";
        }
        cout<<endl;
    }
    return 0;
}