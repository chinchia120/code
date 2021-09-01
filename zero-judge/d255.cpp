#include <iostream>
using namespace std;
long long int GCD(int i,int j){
    while(j%i!=0){
        int tmp=j%i;
        j=i;
        i=tmp;
    }
    return i;
}
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        if(n==0){
            return 0;
        }
        long long int G=0;
        for(int i=1;i<n;i++){
            for(int j=i+1;j<=n;j++){
                G+=GCD(i,j);
            }
        }
        printf("%ld\n",G);
    }
}