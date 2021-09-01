#include <iostream>
using namespace std;
int main(void){
    int year;
    while(cin>>year){
        if(year>0){
            if(year%12==1){
                cout<<"鼠"<<endl;
            }else if(year%12==2){
                cout<<"牛"<<endl;
            }else if(year%12==3){
                cout<<"虎"<<endl;
            }else if(year%12==4){
                cout<<"兔"<<endl;
            }else if(year%12==5){
                cout<<"龍"<<endl;
            }else if(year%12==6){
                cout<<"蛇"<<endl;
            }else if(year%12==7){
                cout<<"馬"<<endl;
            }else if(year%12==8){
                cout<<"羊"<<endl;
            }else if(year%12==9){
                cout<<"猴"<<endl;
            }else if(year%12==10){
                cout<<"雞"<<endl;
            }else if(year%12==11){
                cout<<"狗"<<endl;
            }else if(year%12==0){
                cout<<"豬"<<endl;
            }
        }else if(year<0){
            year=-year;
            if(year%12==1){
                cout<<"豬"<<endl;
            }else if(year%12==2){
                cout<<"狗"<<endl;
            }else if(year%12==3){
                cout<<"雞"<<endl;
            }else if(year%12==4){
                cout<<"猴"<<endl;
            }else if(year%12==5){
                cout<<"羊"<<endl;
            }else if(year%12==6){
                cout<<"馬"<<endl;
            }else if(year%12==7){
                cout<<"蛇"<<endl;
            }else if(year%12==8){
                cout<<"龍"<<endl;
            }else if(year%12==9){
                cout<<"兔"<<endl;
            }else if(year%12==10){
                cout<<"虎"<<endl;
            }else if(year%12==11){
                cout<<"牛"<<endl;
            }else if(year%12==0){
                cout<<"鼠"<<endl;
            }
        }
    }
}