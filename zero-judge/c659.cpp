#include <iostream>
using namespace std;
int main(void){
    string str1,str2;
    cin>>str1;
    getline(cin,str2);
    for(int i=0;i<str2.length();i++){
        if(i==0){
            cout<<str2[i];
        }else{
            if(isspace(str2[i])){
                cout<<" "<<str1<<" ";
            }else{
                cout<<str2[i];
            }
        }
    }
    cout<<endl;
    return 0;
}