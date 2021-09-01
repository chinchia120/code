#include <iostream>
using namespace std;
int main(void){
    int arr[10];
    for(int i=0;i<10;i++){
        cin>>arr[i];
    }
    int hei;
    cin>>hei;
    int cnt=0;
    for(int i=0;i<10;i++){
        if(hei+30>=arr[i]){
            cnt++;
        }
    }
    cout<<cnt<<endl;
}