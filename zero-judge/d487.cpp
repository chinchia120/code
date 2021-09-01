#include <iostream>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        printf("%d! = ",n);
        if(n==0 || n==1){
            printf("1 = 1\n");
            continue;
        }
        long long int sum=1;
        for(int i=n;i>0;i--){
            sum*=i;
            printf("%d ",i);
            if(i!=1){
                printf("* ");
            }
        }
        printf("= %ld\n",sum);
    }
}