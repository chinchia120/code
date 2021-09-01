#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            int a,b;
            cin>>a>>b;
            if(a<b){
                cout<<"impossible\n";
            }else{
                if((a+b)%2!=0 || (a-b)%2!=0){
                    cout<<"impossible\n";
                }else{
                    cout<<(a+b)/2<<" "<<(a-b)/2<<endl;
                }
            }
        }
    }
}