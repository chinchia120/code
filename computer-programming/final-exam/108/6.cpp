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

int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        for(int i=0;i<n;i++){
            char str[4][4];//={{'#','#','#','#'},{'#','2','3','1'},{'#','3','2','4'},{'2','#','#','3'}};
            int cnt=0;
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    cin>>str[j][k];
                    if(str[j][k]=='#'){
                        cnt++;
                    }
                }
            }
            //cout<<cnt<<endl;
            while(cnt!=(0)){
                for(int j=0;j<4;j++){
                    for(int k=0;k<4;k++){
                        if(str[j][k]=='#'){
                            int sumB=0,sumV=0,sumH=0;
                            if(j<2 && k<2){
                                for(int l=0;l<2;l++){
                                    for(int m=0;m<2;m++){
                                        if(str[l][m]!='#'){
                                            sumB=sumB+(str[l][m]-'0');
                                        }
                                    }
                                }
                                for(int l=0;l<4;l++){
                                    if(str[l][k]!='#'){
                                        sumV=sumV+(str[l][k]-'0');
                                    }
                                }
                                for(int l=0;l<4;l++){
                                    if(str[j][l]!='#'){
                                        sumH=sumH+(str[j][l]-'0');
                                    }
                                }
                                for(int l=1;l<5;l++){
                                    if((sumH+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }else if((sumV+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }else if((sumB+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }
                                }
                            }else if(j<2 && k<4){
                                for(int l=0;l<2;l++){
                                    for(int m=2;m<4;m++){
                                        if(str[l][m]!='#'){
                                            sumB=sumB+(str[l][m]-'0');
                                        }
                                    }
                                }
                                for(int l=0;l<4;l++){
                                    if(str[l][k]!='#'){
                                        sumV=sumV+(str[l][k]-'0');
                                    }
                                }
                                for(int l=0;l<4;l++){
                                    if(str[j][l]!='#'){
                                        sumH=sumH+(str[j][l]-'0');
                                    }
                                }
                                for(int l=1;l<5;l++){
                                    if((sumH+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }else if((sumV+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }else if((sumB+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }
                                }
                            }else if(j<4 && k<2){
                                for(int l=2;l<4;l++){
                                    for(int m=0;m<2;m++){
                                        if(str[l][m]!='#'){
                                            sumB=sumB+(str[l][m]-'0');
                                        }
                                    }
                                }
                                for(int l=0;l<4;l++){
                                    if(str[l][k]!='#'){
                                        sumV=sumV+(str[l][k]-'0');
                                    }
                                }
                                for(int l=0;l<4;l++){
                                    if(str[j][l]!='#'){
                                        sumH=sumH+(str[j][l]-'0');
                                    }
                                }
                                for(int l=1;l<5;l++){
                                    if((sumH+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }else if((sumV+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }else if((sumB+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }
                                }
                            }else{
                                for(int l=2;l<4;l++){
                                    for(int m=2;m<4;m++){
                                        if(str[l][m]!='#'){
                                            sumB=sumB+(str[l][m]-'0');
                                        }
                                    }
                                }
                                for(int l=0;l<4;l++){
                                    if(str[l][k]!='#'){
                                        sumV=sumV+(str[l][k]-'0');
                                    }
                                }
                                for(int l=0;l<4;l++){
                                    if(str[j][l]!='#'){
                                        sumH=sumH+(str[j][l]-'0');
                                    }
                                }
                                for(int l=1;l<5;l++){
                                    if((sumH+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }else if((sumV+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }else if((sumB+l)==10){
                                        str[j][k]=l+'0';
                                        cnt--;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            //cout<<cnt<<endl;
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    cout<<str[j][k];
                }
                cout<<endl;
            }
        }
    }
    return 0;
}