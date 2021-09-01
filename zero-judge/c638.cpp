#include <iostream>
using namespace std;
int main(void){
    int year;
    while(cin>>year){
        while(1){
            if(year%10==1){
                cout<<"辛";
                break;
            }else if(year%10==2){
                cout<<"壬";
                break;
            }else if(year%10==3){
                cout<<"癸";
                break;
            }else if(year%10==4){
                cout<<"甲";
                break;
            }else if(year%10==5){
                cout<<"乙";
                break;
            }else if(year%10==6){
                cout<<"丙";
                break;
            }else if(year%10==7){
                cout<<"丁";
                break;
            }else if(year%10==8){
                cout<<"戊";
                break;
            }else if(year%10==9){
                cout<<"己";
                break;
            }else if(year%10==0){
                cout<<"庚";
                break;
            }
        }
        while(1){    
            if(year%12==1){
                cout<<"酉"<<endl;
                break;
            }else if(year%12==2){
                cout<<"戌"<<endl;
                break;
            }else if(year%12==3){
                cout<<"亥"<<endl;
                break;
            }else if(year%12==4){
                cout<<"子"<<endl;
                break;
            }else if(year%12==5){
                cout<<"丑"<<endl;
                break;
            }else if(year%12==6){
                cout<<"寅"<<endl;
                break;
            }else if(year%12==7){
                cout<<"卯"<<endl;
                break;
            }else if(year%12==8){
                cout<<"辰"<<endl;
                break;
            }else if(year%12==9){
                cout<<"巳"<<endl;
                break;
            }else if(year%12==10){
                cout<<"午"<<endl;
                break;
            }else if(year%12==11){
                cout<<"未"<<endl;
                break;
            }else if(year%12==0){
                cout<<"申"<<endl;
                break;
            }
        }    
    }
    return 0;
}