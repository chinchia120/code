#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            string str;
            cin>>str;
            if(str.length()==5){
                cout<<"3\n";
            }else{
                if((str[0]=='o' && str[1]=='n')||(str[0]=='o' && str[2]=='e')||(str[1]=='n' && str[2]=='e')){
                    cout<<"1\n";
                }else{
                    cout<<"2\n";
                }
            }
        }
    }
    return 0;
}