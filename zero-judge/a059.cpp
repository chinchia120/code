#include <iostream>
#include <math.h>
using namespace std;

int main(){
    int n,a,b;
    while(cin>>n){
        for(int i=1;i<n+1;i++){
            cin>>a>>b;
            int sum=0;
            for(int j=a;j<b+1;j++){
                for(int k=1;k<sqrt(j)+1;k++){
                    if(j==k*k){
                        sum=sum+j;
                    }
                }
            }
            cout<<"Case "<<i<<": "<<sum<<endl;
        }
    }
    return 0;
}