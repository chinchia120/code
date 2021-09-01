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
    string ans;
    while(cin>>ans){
        int n;
        cin>>n;
        for(int i=0;i<n;i++){
            string gu;
            cin>>gu;
            int A=0,B=0;
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    if(j==k && ans[j]==gu[k]){
                        A++;
                    }else if(j!=k && ans[j]==gu[k]){
                        B++;
                    }
                }
            }
            printf("%dA%dB\n",A,B);
        }
    }
    return 0;
}