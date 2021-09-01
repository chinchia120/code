#include <iostream>
#include <string>
using namespace std;
int main(void){
    string str;
    getline(cin,str);
    for(int i=0;i<str.length();i++){
        if(i==0){
            cout<<(char)(str[i]-32);
        }else{
            if(isspace(str[i])){
                cout<<endl;
                cout<<(char)(str[i+1]-32);
                i++;
            }else{
                cout<<str[i];
            }
        }
    }
    return 0;
}