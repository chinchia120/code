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
void rf(int a,int b){
    int tmp=b;
    for(int i=2;i<=tmp;i++){
        while(a%i==0 && b%i==0){
            a=a/i;
            b=b/i;
        }
    }
    cout<<a<<"/"<<b<<endl;
}
int main() {
    int n;
    while(scanf("%d",&n)!=EOF){
        for(int i=0;i<n;i++){
            int a,b,c,d;
            cin>>a>>b>>c>>d;
            int nume=a*d+b*c,deno=b*d;
            if(nume/deno<0){
                cout<<"-";
                nume=abs(nume),deno=abs(deno);
            }
            if(nume%deno==0){
                cout<<nume/deno<<endl;
            }else if(nume/deno>=1){
                cout<<nume/deno<<" ";
                nume=nume-(nume/deno)*deno;
                rf(nume,deno);
            }else{
                rf(nume,deno);
            }
        }
    }
    return 0;
}
