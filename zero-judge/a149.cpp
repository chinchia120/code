#include <iostream>
using namespace std;

int main(){
    int n;
    string str;
    while(cin>>n){
        for(int i=0;i<n;i++){
            int sum=1;
            cin>>str;
            for(int j=0;j<str.length();j++){
                sum=sum*(char)(str[j]-48);
            }
            cout<<sum<<endl;
        }
    }
    return 0;
}