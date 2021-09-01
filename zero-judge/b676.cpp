#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        if(n%5==0){
            cout<<"U"<<endl;
        }else if(n%5==1){
            cout<<"G"<<endl;
        }else if(n%5==2){
            cout<<"Y"<<endl;
        }else if(n%5==3){
            cout<<"T"<<endl;
        }else{
            cout<<"I"<<endl;
        }
    }
}