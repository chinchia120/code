#include <iostream>
using namespace std;
int main(void){
    int n,m;
    while(scanf("%d %d",&n,&m)!=EOF){
        int arr[n];
        for(int i=0;i<n;i++){
            scanf("%d",&arr[i]);
        }
        for(int i=0;i<m;i++){
        	int a;
            scanf("%d",&a);
            int flag=0;
            for(int j=0;j<n;j++){
                if(arr[j]==a){
                    printf("%d\n",j+1);
                    flag=1;
                    break;
                }
            }
            if(flag==0){
                printf("0\n");
                continue;
            }
        }
    }
    return 0;
}