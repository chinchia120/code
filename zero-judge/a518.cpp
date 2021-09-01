#include <iostream>
using namespace std;
int main(void){
    int a,b;
    while(cin>>a>>b){
        if(a==(-1) && b==(-1)){
            return 0;
        }else{
            int max,min;
            if(a>b){
                max=a;
                min=b;
            }else{
                max=b;
                min=a;
            }
            if(max-min>50){
                cout<<100-(max-min)<<endl;
            }else{
                cout<<(max-min)<<endl;
            }
        }
    }
}