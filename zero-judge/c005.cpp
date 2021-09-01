#include <iostream>
using namespace std;
int main(void){
    int n;
    scanf("%d",&n);
    for(int i=0;i<n;i++){
        int m;
        cin>>m;
        long long int sum=0;
        for(int j=0;j<m;j++){
            long long int a,b,c;
            cin>>a>>b>>c;
            sum=sum+a*c;
        }
        cout<<sum<<endl;
    }
    return 0;
}