#include <iostream>
#include <string>
using namespace std;
int main(void){
    string str;
    while(cin>>str){
        if(str=="0"){
            return 0;
        }
        int sum=0;
        for(int i=0;i<str.length();i++){
            if(i%2==0){
                sum+=str[i]-'0';
            }else{
                sum-=str[i]-'0';
            }
        }
        if(sum%11==0){
            cout<<str<<" is a multiple of 11."<<endl; 
        }else{
            cout<<str<<" is not a multiple of 11."<<endl;
        }
    }
    return 0;
}