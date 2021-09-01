#include <iostream>
using namespace std;

int main(){
    int a,b,c,sum=0;
    for(int i=0;i<5;i++){
        cin>>a>>b>>c;
        if(a+b>c&&a+c>b&&b+c>a){
            sum=sum+1;
        }
    }
    cout<<sum<<endl;
    return 0;
}