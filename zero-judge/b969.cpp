#include <iostream>
#include <cstring>
using namespace std;
int main(void){
    string str1,str2;
    getline(cin,str1);
    getline(cin,str2);
    for(int i=0;i<str1.length();i++){
        if(i==0){
            cout<<str2<<", "<<str1[i];
        }else{
            if(isspace(str1[i])){
                cout<<endl;
                cout<<str2<<", ";
            }else{
                cout<<str1[i];
            }
        }
    }
    return 0;
}