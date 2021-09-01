#include <iostream>
#include <string>
using namespace std;
int main(void){
    string str1,str2;
    int a=0,b=0;
    while(cin>>str1>>str2){
        if((str1=="Scissors" && str2=="Stone") || (str1=="Stone" && str2=="Paper") || (str1=="Paper" && str2=="Scissors")){
            cout<<"靈夢獲勝"<<endl;
            a++;
        }else if((str2=="Scissors" && str1=="Stone") || (str2=="Stone" && str1=="Paper") || (str2=="Paper" && str1=="Scissors")){
            cout<<"紫獲勝"<<endl;
            b++;
        }else{
            if(a>b){
                cout<<"悲慘的籌措起香油錢"<<endl;
            }else{
                cout<<"螢火的蹤跡"<<endl;
            }
        }
    }
}