#include <iostream>
#include <string>
using namespace std;
int main(void){
    string a,b;
    while(cin>>a>>b){
        if(a=="0" && b=="0"){
            return 0;
        }
        int len;
        if(a.length()>b.length()){
            len=a.length();
        }else{
            len=b.length();
        }
        int arrA[len]={0},arrB[len]={0};
        for(int i=0;i<a.length();i++){
            arrA[i]=a.at(a.length()-1-i)-'0';
        }
        for(int i=0;i<b.length();i++){
            arrB[i]=b.at(b.length()-1-i)-'0';
        }
        int cnt=0,tmp;
        for(int i=0;i<len;i++){
            if(i==0){
                if(arrA[i]+arrB[i]>=10){
                    tmp=(arrA[i]+arrB[i])/10;
                    cnt++;
                }else{
                    tmp=0;
                }
            }else{
                if(arrA[i]+arrB[i]+tmp>=10){
                    tmp=(arrA[i]+arrB[i]+tmp)/10;
                    cnt++;
                }else{
                    tmp=0;
                }
            }
        }
        if(cnt==0){
            printf("No carry operation.\n");
        }else if(cnt==1){
            printf("1 carry operation.\n");
        }else{
            printf("%d carry operations.\n",cnt);
        }
    }
}