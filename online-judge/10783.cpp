#include <iostream>
using namespace std;
int main(void) {
    int n;
    scanf("%d",&n);
    for(int i=0;i<n;i++){
        int start,end;
        cin>>start>>end;
        if(start%2==0){
            start=start+1;
        }
        if(end%2==0){
            end=end-1;
        }
        long long int sum=0;
        for(int j=start;j<end+1;j+=2){
            sum=sum+j;
        }
        cout<<"Case "<<i+1<<": "<<sum<<endl;
    }
    return 0;
}