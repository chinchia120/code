#include <iostream>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        int max=0;
        for(int i=0;i<n;i++){
            int num;
            cin>>num;
            if(num>max){
                max=num;
            }
        }
        cout<<max<<endl;
    }
}