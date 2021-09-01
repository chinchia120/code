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
            int sum=0;
            for(int j=0;j<str.length();j++){
                if(str[j]=='a' || str[j]=='e' || str[j]=='i' || str[j]=='o' || str[j]=='u'){
                    sum++;
                    if(str[j+1]=='a' || str[j+1]=='e' || str[j+1]=='i' || str[j+1]=='o' || str[j+1]=='u'){
                        sum--;
                    }
                }
            }
            cout<<sum<<endl;
        }
    }
    return 0;
}