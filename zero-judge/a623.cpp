#include <iostream>
using namespace std;
int count_nume(int n,int m){
    int cnt=0,sum=1;
    while(cnt!=m){
        sum=sum*n;
        cnt++;
        n--;
    }
    return sum;
}
int count_deno(int m){
    int sum=1;
    for(int i=1;i<=m;i++){
        sum=sum*i;
    }
    return sum;
}
int main(void){
    int n,m;
    while(cin>>n>>m){
        int nume=count_nume(n,m);
        int deno=count_deno(m);
        //cout<<nume<<" "<<deno<<endl;
        cout<<nume/deno<<endl;
    }
    return 0;
}