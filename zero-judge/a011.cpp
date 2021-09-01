#include <iostream>
#include <string>
using namespace std;
int main(void){
    string str;
    while(getline(cin,str)){
        int cnt=0;
        for(int i=0;i<str.length()-1;i++){
            if(isalpha(str.at(i))){
                cnt++;
                if(isalpha(str.at(i+1))){
                    cnt--;
                }
            }
        }
        cout<<cnt<<endl;
    }
    return 0;
}