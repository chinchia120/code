#include <iostream>
using namespace std;

int main(){
    int a,b,sum=0;
    while(cin>>a>>b){
        for(int i=a;i<b+1;i++){
            if(i%2==0){
                sum=sum+i;
            }
        }
        cout<<sum<<endl;
    }
    return 0;
}