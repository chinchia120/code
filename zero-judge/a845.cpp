#include <iostream>
using namespace std;
int main(void){
    int n;
    cin>>n;
    int num[n];
    for(int i=0;i<n;i++){
        cin>>num[i];
    }
    int q;
    cin>>q;
    for(int i=0;i<q;i++){
        int n1,n2;
        cin>>n1>>n2;
        cout<<num[n1]+num[n2]<<endl;
    }
}