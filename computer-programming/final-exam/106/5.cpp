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
            string str1,str2;
            cin>>str1>>str2;
            int arr[str2.length()],check=0;
            for(int j=0;j<str1.length();j++){
                for(int k=0;k<str2.length();k++){
                    if(str1[j]==str2[k]){
                        int a=j,b=0,flag=1;
                        while(b<str2.length()){
                            if(str1[a]==str2[b]){
                                a++;
                                b++;
                            }else{
                                flag=0;
                                break;
                            }
                        }
                        if(flag==1){
                            check=1;
                        }
                    }
                }   
            }
            if(check==1){
                cout<<"Yes\n";
            }else{
                cout<<"No\n";
            }
        }
    }
    return 0;
}