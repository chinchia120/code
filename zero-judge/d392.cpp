#include <iostream>
#include <string>
using namespace std;
int main(void){
    string str;
    while(getline(cin,str)){
        long long unsigned int sum=0;
        for(int i=0;i<str.length();i++){
            if(isdigit(str[i])){
                sum=sum+(str[i]-'0');
            }    
        }
        cout<<sum<<endl;
    }
    return 0;
}