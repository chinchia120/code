#include <iostream>
using namespace std;

int main(void){
    int n;
    while(scanf("%d",&n)){
        if(n==0){
            return 0;
        }
        long long int arr[n],sum=1;
        for(int i=0;i<n;i++){
            cin>>arr[i];
            //cout<<arr[i]<<" ";
        }
        for(int i=0;i<n;i++){
            sum=sum*arr[i];
        }
        //cout<<sum<<endl;
        for(long long int i=arr[0];i<sum+1;i++){
            int flag=0;
            for(int j=0;j<n;j++){
                if(i%arr[j]!=0){
                    flag=1;
                }
            }
            if(flag==0){
                printf("%d\n",i);
                break;
            }
        }
    }
}