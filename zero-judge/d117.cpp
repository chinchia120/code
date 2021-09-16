#include <iostream>
#include <math.h>
#include <string>
#include <ctype.h>
using namespace std;
int main(void){
    string str;
    while(cin>>str){
        if(str=="a"){
            printf("It is a prime word.\n");
            continue;
        }
        long long int sum=0;
        for(int i=0;i<str.length();i++){
            if(islower(str[i])){
                sum+=(char)str[i]-96;
            }else{
                sum+=(char)str[i]-38;
            }
        }
        int flag=1;
        for(int i=2;i<sum;i++){
            if(sum%i==0){
                flag=0;
                break;
            }
        }
        if(flag==1){
            printf("It is a prime word.\n");
        }else{
            printf("It is not a prime word.\n");
        }
    }
    return 0;
}