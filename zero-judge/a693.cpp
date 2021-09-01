#include <iostream>
using namespace std;
int main(void){
    int n,m;
    while(scanf("%d %d",&n,&m)!=EOF){
        int arr[n]={};
        for(int i=0;i<n;i++){
            scanf("%d",&arr[i]);
        }
        for(int i=0;i<m;i++){
            int start,end;
            scanf("%d %d",&start,&end);
            int sum=0;
            for(int j=start-1;j<end;j++){
                sum=sum+arr[j];
            }
            printf("%d\n",sum);
        }
    }    
    return 0;
}