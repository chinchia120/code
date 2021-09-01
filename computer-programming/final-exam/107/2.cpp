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
            int a;
            string str;
            cin>>a>>str;
            if(str=="Upper"){
                cout<<(char)(a+64)<<endl;
            }else{
                cout<<(char)(a+96)<<endl;
            }
        }
    }
    return 0;
}