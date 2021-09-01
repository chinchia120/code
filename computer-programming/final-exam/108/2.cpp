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
    for(int i=2;i<tmp;i++){
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
            int a,b;
            cin>>a>>b;
            if(a%b==0){
                cout<<a/b<<endl;
            }else if(a/b>=1){
                cout<<a/b<<" ";
                a=a-(a/b)*b;
                rf(a,b);
            }else{
                rf(a,b);
            }
        }
    }
    return 0;
}