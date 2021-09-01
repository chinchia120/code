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
    string ans;
    while(cin>>ans){
        int n;
        cin>>n;
        for(int i=0;i<n;i++){
            string guess;
            cin>>guess;
            int n=0,m=0;
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    if(j==k && ans[j]==guess[k]){
                        n++;
                    }else if(j!=k && ans[j]==guess[k]){
                        m++;
                    }
                }
            }
            cout<<100*n+10*m<<endl;
        }
    }
    return 0;
}