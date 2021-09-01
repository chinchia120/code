#include <iostream>
#include <vector>
using namespace std;
int main(void){
    string _name;
    vector<string> name;
    while(getline(cin,_name)){
        int flag=0;
        for(int i=0;i<name.size();i++){
            if(_name==name[i]){
                flag=1;
                break;
            }
        }
        if(flag==1){
            cout<<"YES"<<endl;
        }else{
            cout<<"NO"<<endl;
            name.push_back(_name);
        }
    }   
    return 0;
}