#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            double a,b,c;
            cin>>a>>b;
            c=(b-a)/a*100;
            if((100*(b-a)>=10*a) || (100*(b-a)<=-7*a)){
                printf("%04.2f%% dispose\n",c);
            }else{
                printf("%04.2f%% keep\n",c);
            }
        }
    }
    return 0;
}