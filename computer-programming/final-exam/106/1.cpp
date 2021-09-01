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
            int l,w,h;
            cin>>l>>w>>h;
            if(l*w>l*h && l*w>w*h){
                cout<<l*w<<endl;
            }else if(l*h>l*w && l*h>w*h){
                cout<<l*h<<endl;
            }else{
                cout<<w*h<<endl;
            }
        }
    }
    return 0;
}
