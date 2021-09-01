#include <iostream>
using namespace std;
int main(void){
    string str;
    while(getline(cin,str)){
        for(int i=0;i<str.length();i++){
            if(isspace(str[i])){
                cout<<" ";
            }else{
                if(str[i]=='1'){
                    cout<<"`";
                }else if(str[i]=='2'){
                    cout<<"1";
                }else if(str[i]=='3'){
                    cout<<"2";
                }else if(str[i]=='4'){
                    cout<<"3";
                }else if(str[i]=='5'){
                    cout<<"4";
                }else if(str[i]=='6'){
                    cout<<"5";
                }else if(str[i]=='7'){
                    cout<<"6";
                }else if(str[i]=='8'){
                    cout<<"7";
                }else if(str[i]=='9'){
                    cout<<"8";
                }else if(str[i]=='0'){
                    cout<<"9";
                }else if(str[i]=='-'){
                    cout<<"0";
                }else if(str[i]=='='){
                    cout<<"-";
                }else if(str[i]=='W'){
                    cout<<"Q";
                }else if(str[i]=='E'){
                    cout<<"W";
                }else if(str[i]=='R'){
                    cout<<"E";
                }else if(str[i]=='T'){
                    cout<<"R";
                }else if(str[i]=='Y'){
                    cout<<"T";
                }else if(str[i]=='U'){
                    cout<<"Y";
                }else if(str[i]=='I'){
                    cout<<"U";
                }else if(str[i]=='O'){
                    cout<<"I";
                }else if(str[i]=='P'){
                    cout<<"O";
                }else if(str[i]=='['){
                    cout<<"P";
                }else if(str[i]==']'){
                    cout<<"[";
                }else if(str[i]=='\\'){
                    cout<<"]";
                }else if(str[i]=='S'){
                    cout<<"A";
                }else if(str[i]=='D'){
                    cout<<"S";
                }else if(str[i]=='F'){
                    cout<<"D";
                }else if(str[i]=='G'){
                    cout<<"F";
                }else if(str[i]=='H'){
                    cout<<"G";
                }else if(str[i]=='J'){
                    cout<<"H";
                }else if(str[i]=='K'){
                    cout<<"J";
                }else if(str[i]=='L'){
                    cout<<"K";
                }else if(str[i]==';'){
                    cout<<"L";
                }else if(str[i]=='\''){
                    cout<<";";
                }else if(str[i]=='X'){
                    cout<<"Z";
                }else if(str[i]=='C'){
                    cout<<"X";
                }else if(str[i]=='V'){
                    cout<<"C";
                }else if(str[i]=='B'){
                    cout<<"V";
                }else if(str[i]=='N'){
                    cout<<"B";
                }else if(str[i]=='M'){
                    cout<<"N";
                }else if(str[i]==','){
                    cout<<"M";
                }else if(str[i]=='.'){
                    cout<<",";
                }else if(str[i]=='/'){
                    cout<<".";
                }       
            }
        }
        cout<<endl;
    }
    return 0;
}