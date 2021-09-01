#include <iostream>
#include <string>
#include <cmath>
using namespace std;
int transform(int n,string str){
    int k=1,sum=0;
    for(int i=str.length()-1;i>=0;i--){
        sum=sum+(str[i]-'0')*k;
        k=k*n;
    }
    return sum;
}
int digit(int ans){
    int cnt=0;
    while(ans!=0){
        cnt++;
        ans=ans/10;
    }
    return cnt;
}
int main(void){
    int n;
    string str;
    cin>>n>>str;
    int ans=transform(n,str);
    int dig=str.length();
    int sum=0;
    for(int i=0;i<str.length();i++){
        sum=sum+pow(str[i]-'0',dig);
    }
    if(sum==ans){
        printf("YES\n");
    }else{
        printf("NO\n");
    }
    return 0;
}