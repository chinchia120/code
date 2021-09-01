#include <iostream>
using namespace std;
int main(void){
    int occupation,level;
    while(scanf("%d %d",&occupation,&level)!=EOF){
        if(occupation==0){
            cout<<"0"<<endl;
        }else if(occupation==1 || occupation==3 || occupation==4){
            if(level<10){
                cout<<"0"<<endl;
            }else if(level<30){
                cout<<3*(level-10)+1<<endl;
            }else if(level<70){
                cout<<3*(level-10)+2<<endl;
            }else if(level<120){
                cout<<3*(level-10)+3<<endl;
            }else{
                cout<<3*(level-10)+6<<endl;
            }
        }else{
            if(level<8){
                cout<<"0"<<endl;
            }else if(level<30){
                cout<<3*(level-8)+1<<endl;
            }else if(level<70){
                cout<<3*(level-8)+2<<endl;
            }else if(level<120){
                cout<<3*(level-8)+3<<endl;
            }else{
                cout<<3*(level-8)+6<<endl;
            }
        }
    }
    return 0;
}