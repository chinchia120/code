#include <iostream>
using namespace std;

int main(){
    long long int n,m,i=1;
    while(cin>>n>>m){
        if(m==0){
            cout<<"Go Kevin!!"<<endl;
            continue;
        }
        for(int i=1;i<9999999999;i=i+m){
            n=n-i;
            if(n==0){
                cout<<"Go Kevin!!"<<endl;
                break;
            }else if(n<0){
                cout<<"No Stop!!"<<endl;
                break;
            }
        }
    }
    return 0;
}