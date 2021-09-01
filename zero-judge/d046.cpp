#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        int cnt=0;
        for(int i=0;i<n;i++){
            int weight;
            cin>>weight;
            if(weight<=10){
                cnt++;
            }
        }
        cout<<cnt<<endl;
    }
    return 0;
}