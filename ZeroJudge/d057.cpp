#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
using namespace std;
int main(void){
    int x1,y1,x2,y2;
    while(cin>>x1>>y1>>x2>>y2){
        if(x1==0 && y1==0 && x2==0 && y2==0){
            return 0;
        }else if(x1==x2 && y1==y2){
            cout<<"0"<<endl;
        }else if(x1==x2 && y1!=y2){
            cout<<"1"<<endl;
        }else if(x1!=x2 && y1==y2){
            cout<<"1"<<endl;
        }else{
            float p=fabs(x1-x2);
            float q=fabs(y1-y2);
            if(p==q){
                cout<<"1"<<endl;
            }else{
                cout<<"2"<<endl;
            }
        }
    }
}