#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        int tmp=0,kill=0,assist=0,die=0;
        for(int i=0;i<n;i++){
            string str;
            cin>>str;
            if(str=="Get_Kill"){
                tmp++;
                kill++;
                if(tmp<3){
                    cout<<"You have slain an enemie."<<endl;
                }else if(tmp==3){
                    cout<<"KILLING SPREE!"<<endl;
                }else if(tmp==4){
                    cout<<"RAMPAGE~"<<endl;
                }else if(tmp==5){
                    cout<<"UNSTOPPABLE!"<<endl;
                }else if(tmp==6){
                    cout<<"DOMINATING!"<<endl;
                }else if(tmp==7){
                    cout<<"GUALIKE!"<<endl;
                }else if(tmp>7){
                    cout<<"LEGENDARY!"<<endl;
                }
            }else if(str=="Get_Assist"){
                assist++;
            }else if(str=="Die"){
                die++;
                if(tmp<3){
                    cout<<"You have been slained."<<endl;
                    tmp=0;
                }else{
                    cout<<"SHUTDOWN."<<endl;
                    tmp=0;
                }
            }
        }
        cout<<kill<<"/"<<die<<"/"<<assist<<endl;
    }
    return 0;
}