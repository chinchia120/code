#include <iostream>
using namespace std;
int main(void){
    int n,cnt=1;
    while(cin>>n){
        if(n==0){
            return 0;
        }
        int arr[n];
        long long int sum=0,average;
        for(int i=0;i<n;i++){
            cin>>arr[i];
            sum=sum+arr[i];
        }
        average=sum/n;
        sum=0;
        for(int i=0;i<n;i++){
            if(arr[i]<average){
                sum=sum+average-arr[i];
            }
        }
        printf("Set #%d\n",cnt);
        printf("The minimum number of moves is %ld.\n\n",sum);
        cnt++;
    }
}