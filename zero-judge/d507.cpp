#include <iostream>
using namespace std;

int main(){
    float a,b,c,cos,l,m,s;
    while(cin>>a>>b>>c){
        if(a>b&&a>c){
            if(b>c||b==c){
                l=a;
                m=b;
                s=c;
            }else{
                l=a;
                m=c;
                s=b;
            }   
        }else if(b>a&&b>c){
            if(a>c||a==c){
                l=b;
                m=a;
                s=c;
            }else{
                l=b;
                m=c;
                s=a;
            }
        }else if(c>a&&c>b){
            if(a>b||a==b){
                l=c;
                m=a;
                s=b;
            }else{
                l=c;
                m=b;
                s=a;
            }
        }
        cos=(s*s+m*m-l*l)/(2*s*m);
        if(cos>0){
            cout<<"acute triangle"<<endl;
        }else if(cos==0){
            cout<<"right triangle"<<endl;
        }else if(cos<0){
            cout<<"obtuse triangle"<<endl;
        }
    }
    return 0;
}