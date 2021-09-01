#include <iostream>
using namespace std;

int main(){
    int m,n,sum,a,b;
    while(cin>>m>>n){
        sum=m*60+n+150;
        if(sum>1439){
            a=(sum-1440)/60;
            b=(sum-1440)%60;
            if(a<10){
                if(b==0){
                    cout<<"0"<<a<<":"<<"00"<<endl;
                }else if(b<10){
                    cout<<"0"<<a<<":"<<"0"<<b<<endl;
                }else{
                    cout<<"0"<<a<<":"<<b<<endl;
                }
            }else{
                if(b==0){
                    cout<<a<<":"<<"00"<<endl;
                }else if(b<10){
                    cout<<a<<":"<<"0"<<b<<endl;
                }else{
                    cout<<a<<":"<<b<<endl;
                }
            }
        }else{
            a=sum/60;
            b=sum%60;
            if(a<10){
                if(b==0){
                    cout<<"0"<<a<<":"<<"00"<<endl;
                }else if(b<10){
                    cout<<"0"<<a<<":"<<"0"<<b<<endl;
                }else{
                    cout<<"0"<<a<<":"<<b<<endl;
                }
            }else{
                if(b==0){
                    cout<<a<<":"<<"00"<<endl;
                }else if(b<10){
                    cout<<a<<":"<<"0"<<b<<endl;
                }else{
                    cout<<a<<":"<<b<<endl;
                }
            }
        }
    }
    return 0;
}