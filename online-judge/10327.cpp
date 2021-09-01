#include <iostream>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        int arr[n];
        for(int i=0;i<n;i++){
            cin>>arr[i];
        }
        int cnt=0;
        for(int i=0;i<n;i++){
            for(int j=i+1;j<n;j++){
                if(arr[i]>arr[j]){
                    int tmp=arr[i];
                    arr[i]=arr[j];
                    arr[j]=tmp;
                    cnt++;
                }
            }
        }
        printf("Minimum exchange operations : %d\n",cnt);
    }
    return 0;
}