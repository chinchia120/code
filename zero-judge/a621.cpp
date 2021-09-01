#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        int sum=1;
        for(int i=0;i<n+1;i++){
            printf("2^%d = %d\n",i,sum);
            sum*=2;
        }
    }
}