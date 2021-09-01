#include <iostream>
using namespace std;
int main(void){
    string str;
    long long unsigned int n;
    while(cin>>str>>n){
        long long unsigned int a=0,sum=0;
        if(str=="Sunday"){
            a=1;
        }else if(str=="Monday"){
            a=2;
        }else if(str=="Tuesday"){
            a=3;
        }else if(str=="Wednesday"){
            a=4;
        }else if(str=="Thursday"){
            a=5;
        }else if(str=="Friday"){
            a=6;
        }else if(str=="Saturday"){
            a=7;
        }
        sum=a+n;
        if(sum>=7){
            sum=sum%7;
        }
        if(sum==1){
            cout<<"Sunday"<<endl;
        }else if(sum==2){
            cout<<"Monday"<<endl;
        }else if(sum==3){
            cout<<"Tuesday"<<endl;
        }else if(sum==4){
            cout<<"Wednesday"<<endl;
        }else if(sum==5){
            cout<<"Thursday"<<endl;
        }else if(sum==6){
            cout<<"Friday"<<endl;
        }else if(sum==0){
            cout<<"Saturday"<<endl;
        }
    }
    return 0;
}