#include <iostream>
#include <string>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        string str;
        cin>>str;
        for(int i=0;i<str.length();i++){
            if(islower(str[i])){
                cout<<(char)(str[i]-32);
            }else{
                cout<<(char)(str[i]+32);
            }
        }
        cout<<endl;
    }
}