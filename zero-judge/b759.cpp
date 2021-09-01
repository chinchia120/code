#include <iostream>
using namespace std;
int main(void){
    string str;
    while(cin>>str){
        for(long long int i=0;i<str.length();i++){
            for(long long int j=0;j<str.length();j++){
                cout<<str[(i+j)%str.length()];
            }
            cout<<endl;
        }

    }
    return 0;
}