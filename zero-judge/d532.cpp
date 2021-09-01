#include <iostream>
using namespace std;

int main(){
    int a,b,sum=0;
    while(cin>>a>>b){
        for(int i=a;i<b+1;i++){
            if((i%400==0)||((i%4==0)&&(i%100!=0))){
                sum=sum+1;
            }
        }
        cout<<sum<<endl;
    }
    return 0;
}