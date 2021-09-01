#include <iostream>
#include <string>
using namespace std;
int main(void){
    string str;
    while(cin>>str){
        if(str=="0"){
            return 0;
        }
        int sum=0,flag=1;
        for(int i=0;i<str.length();i++){
            if(isupper(str[i])){
                sum=sum+(str.at(i)-64);
            }else if(islower(str[i])){
                sum=sum+(str.at(i)-96);
            }else{
                cout<<"Fail"<<endl;
                flag=0;
                break;
            }
        }
        if(flag==1){
            cout<<sum<<endl;
        }
    }
}