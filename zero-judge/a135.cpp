#include <iostream>
using namespace std;
int main(void){
    string str;
    int cnt=1;
    while(cin>>str){
        if(str=="#"){
            return 0;
        }
        if(str=="HELLO"){
            printf("Case %d: ENGLISH\n",cnt);
        }else if(str=="HOLA"){
            printf("Case %d: SPANISH\n",cnt);
        }else if(str=="HALLO"){
            printf("Case %d: GERMAN\n",cnt);
        }else if(str=="BONJOUR"){
            printf("Case %d: FRENCH\n",cnt);
        }else if(str=="CIAO"){
            printf("Case %d: ITALIAN\n",cnt);
        }else if(str=="ZDRAVSTVUJTE"){
            printf("Case %d: RUSSIAN\n",cnt);
        }else{
            printf("Case %d: UNKNOWN\n",cnt);
        }
        cnt++;
    }
}