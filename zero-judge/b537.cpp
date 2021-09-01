#include <iostream>
using namespace std;
double F(double a){
    for(int i=i;i<99999;i++){
        if(i==1){
            if(a==i){
                return i;
            }
        }else if(i%2==1){
            if(a==(1/F(a-1))){
                return i;
            }
        }else{
            if(a==(1+F(a/2))){
                return i;
            }
        }
    }
}
int main(void){
    double a,b,c;
    while(cin>>a>>b){
        c=a/b;
        cout<<F(c)<<endl;
    }
    return 0;
}