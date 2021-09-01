#include <iostream>
using namespace std;
int F(int f){
    if(f==1){
        return 1;
    }else if(f%2==0){
        return F(f/2);
    }else{
        return F(f-1)+F(f+1);
    }
}
int main(void){
    int n;
    while(cin>>n){
        cout<<F(n)<<endl;
    }
    return 0;
}