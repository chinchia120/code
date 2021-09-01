#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <stack>
#include <algorithm>
#include <vector>
using namespace std;

int main() {
    int n;
    while(scanf("%d",&n)!=EOF){
        for(int i=0;i<n;i++){
            int med;
            string str;
            cin>>med>>str;
            int A=0,B=0,C=0;
            for(int j=0;j<str.length();j++){
                if(str[j]=='1' || str[j]=='2' || str[j]=='3'){
                    str[j]='a';
                    A++;
                }else if(str[j]=='4' || str[j]=='5'){
                    str[j]='b';
                    B++;
                }else{
                    str[j]='c';
                    C++;
                }
            }
            if(med<=A){
                for(int j=0;j<str.length();j++){
                    if(med==0){
                        cout<<"N";
                        continue;
                    }
                    if(str[j]=='a'){
                        cout<<"Y";
                        med--;
                    }else{
                        cout<<"N";
                    }
                }
            }else if(med<=(A+B)){
                int medB=med-A;
                for(int j=0;j<str.length();j++){
                    if(med==0){
                        cout<<"N";
                        continue;
                    }
                    if(str[j]=='a'){
                        cout<<"Y";
                        med--;
                    }else if(str[j]=='b'){
                        if(medB>0){
                            cout<<"Y";
                            medB--;
                            med--;
                        }else{
                            cout<<"N";
                        }
                    }else{
                        cout<<"N";
                    }
                }
            }else{
                int medC=med-A-B;
                for(int j=0;j<str.length();j++){
                    if(med==0){
                        cout<<"N";
                        continue;
                    }
                    if(str[j]=='a' || str[j]=='b'){
                        cout<<"Y";
                        med--;
                    }else if(str[j]=='c'){ 
                        if(medC>0){
                            cout<<"Y";
                            medC--;
                            med--;
                        }else{
                            cout<<"N";
                        }     
                    }
                }
            }
            cout<<endl;
        }
    }
    return 0;
}
