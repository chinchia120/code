#include <iostream>
using namespace std;

int main(){
    long long int d,sum=0;
    while(cin>>d){
        for(int i=1;i<51;i++){
            sum=sum+(i+d*(i-1)*(i-2)/2);
        }
        cout<<sum<<endl;
    }
    return 0;
}