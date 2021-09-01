#include <iostream>
using namespace std;
int F(int a){
    if(a==1){
        return 1;
    }else if(a==2){
        return 2;
    }else{
        return F(a-1)+F(a-2);
    }
}
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        printf("%d\n",F(n));
    }
    return 0;
}