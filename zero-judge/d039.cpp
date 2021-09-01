#include <iostream>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        for(int i=0;i<n;i++){
            int l,w;
            cin>>l>>w;
            cout<<(l/3)*(w/3)<<endl;
        }
    }
    return 0;
}