#include <iostream>
using namespace std;
int main(void){
    long long int hr,min,sec;
    while(cin>>hr>>min>>sec){
        long long int h=0,m=0,s=0;
        while(sec>=60){
            s=sec/60;
            sec=sec%60;
            min=min+s;
            break;    
        }
        while(min>=60){
            m=min/60;
            min=min%60;
            hr=hr+m;
            break;
        }
        while(hr>=24){
            hr=hr%24;
            break;
        }
        printf("%02d:%02d:%02d\n",hr,min,sec);
        
    }
    return 0;
}