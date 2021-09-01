#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        int x=0,y=0;
        for(int i=0;i<n;i++){
            int a,b;
            cin>>a>>b;
            //cout<<x<<" "<<y<<endl;
            if(a==0){
                y=y+b;
            }else if(a==1){
                x=x+b;
            }else if(a==2){
                y=y-b;
            }else if(a==3){
                x=x-b;
            }
        }
        cout<<x<<" "<<y<<endl;
    }
    return 0;
}