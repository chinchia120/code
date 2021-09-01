#include <iostream>
using namespace std;

int main(){
    int a,b,big,small,sum;
    while(cin>>a>>b){
        if(a>b){
            big=a;
            small=b;
        }else{
            big=b;
            small=a;
        }
        for(int i=small;i<big+1;i++){
            if(i%2==0){
                sum=sum+i;
            }
        }
        cout<<sum<<endl;
    }
    return 0;
}