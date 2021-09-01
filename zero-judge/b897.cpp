#include <iostream>
using namespace std;

int main(){
    int x,y;
    while(cin>>x>>y){
        int sum=1;
        for(int i=x;i>y+1;i++){
            sum=sum*i;
            cout<<sum<<endl;
        }
    }
    return 0;
}