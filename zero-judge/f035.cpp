#include <iostream>
using namespace std;
#define dig str.length()
int main(void){
    string str;
    while(getline(cin,str)){
        for(int i=0;i<dig;i++){
            cout<<int(str[i])<<"_"[i==dig-1];
        }
        cout<<"\n";
    }
    return 0;
}