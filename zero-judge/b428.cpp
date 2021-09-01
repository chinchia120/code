#include <iostream>
#include <string>
using namespace std;
int main(void){
    string s1,s2;
    while(cin>>s1>>s2){
        int a=s1.at(0),b=s2.at(0);
        if(a>b){
            cout<<26-(a-b)<<endl;
        }else{
            cout<<(b-a)<<endl;
        }
    }
}