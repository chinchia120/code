#include <iostream>
using namespace std;
int main(void){
    long long int a;
    while(cin>>a){
        int arr[99999],cnt=0;
        while(a!=1){
            if(a%2==0){
                arr[cnt]=0;
                a=a/2;
                cnt++;
            }else if(a%2==1){
                arr[cnt]=1;
                a=a/2;
                cnt++;
            }
        }
        while(a==1){
            cout<<"1";
            for(int i=cnt-1;i>=0;i=i-1){
                cout<<arr[i];
            }
            cout<<endl;
            break;
        }
    }
    return 0;
}