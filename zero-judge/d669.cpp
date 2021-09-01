#include <iostream>
using namespace std;
int main(void){
    int h1,m1,h2,m2;
    while(cin>>h1>>m1>>h2>>m2){
        int time1=h1*60+m1;
        int time2=h2*60+m2;
        if(time1+time2==0){
            return 0;
        }else if(h1==h2){
            if(time1>time2){
                printf("%d\n",24*60-(time1-time2));
            }else{
                printf("%d\n",(time2-time1));
            }
        }else if(h1>h2){
            printf("%d\n",24*60-(time1-time2));
            
        }else if(h1<h2){
            printf("%d\n",(time2-time1));
        }
    }
}