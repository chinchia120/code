#include <iostream>
using namespace std;

int main(){
    int a,b,c;
    while(cin>>a>>b){
        if(a==0 && b==0){
            return 0;
        }
        if(a>b){
            while(a%b!=0){
                c=a%b;
                a=b;
                b=c;
            }
            while(a%b==0){
                cout<<b<<endl;
                break;
            }
        }else if(a==b){
            cout<<a<<endl;
        }else if(a<b){
            while(b%a!=0){
                c=b%a;
                b=a;
                a=c;
            }    
            while(b%a==0){
                cout<<a<<endl;
                break;
            }
        }
    }
    return 0;
}