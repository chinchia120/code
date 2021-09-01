#include <iostream>
using namespace std;
int main(void){
    long long unsigned int a1,q1,a2,q2,a3,q3;
    while(cin>>a1>>q1>>a2>>q2>>a3>>q3){
        for(long long unsigned int i=a1;i<a1*a2*a3;i++){
            if(i%a1==q1 && i%a2==q2 && i%a3==q3){
                cout<<i<<endl;
                break;
            }
        }
    }
    return 0;
}