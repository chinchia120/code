#include <bits\stdc++.h>
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
            int x,y;
            cin>>x>>y;
            string str;
            cin>>str;
            for(int i=0;i<str.length();i++){
                if(i%2==0){
                    if(str[i]=='N'){
                        y=y+(str[i+1]-'0');
                    }else if(str[i]=='S'){
                        y=y-(str[i+1]-'0');
                    }else if(str[i]=='E'){
                        x=x+(str[i+1]-'0');
                    }else{
                        x=x-(str[i+1]-'0');
                    }
                }
            }
            cout<<"("<<x<<","<<y<<")"<<endl;
        }
    }
    return 0;
}
