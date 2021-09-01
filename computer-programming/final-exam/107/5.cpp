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
            string str;
            cin>>str;
            int x=0,y=0;
            string direction="+x";
            for(int j=0;j<str.length();j++){
                if(str[j]=='F'){
                    if(direction=="+x"){
                        x++;
                    }else if(direction=="-x"){
                        x--;
                    }else if(direction=="+y"){
                        y++;
                    }else{
                        y--;
                    }
                }else if(str[j]=='R'){
                    if(direction=="+x"){
                        direction="-y";
                    }else if(direction=="-x"){
                        direction="+y";
                    }else if(direction=="+y"){
                        direction="+x";
                    }else{
                        direction="-x";
                    }
                }else{
                    if(direction=="+x"){
                        direction="+y";
                    }else if(direction=="-x"){
                        direction="-y";
                    }else if(direction=="+y"){
                        direction="-x";
                    }else{
                        direction="+x";
                    }
                }
            }
            printf("(%d,%d)\n",x,y);
        }
    }
    return 0;
}