#include <iostream>
using namespace std;
#define dig str.length()
int main(void){
    string str;
    while(cin>>str){
        cout<<str[0];
        for(int i=0;i<dig-2;i++){
            cout<<"_";
        }
        cout<<str[dig-1]<<endl;
    }
    return 0;
}