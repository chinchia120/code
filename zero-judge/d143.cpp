#include <iostream>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        for(int i=0;i<n;i++){
            long long int a,b;
            scanf("%ld %ld",&a,&b);
            if(a>b){
                printf(">\n");
            }else if(a<b){
                printf("<\n");
            }else{
                printf("=\n");
            }
        }
    }
    return 0;
}