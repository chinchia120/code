#include <iostream>
#include <string>
using namespace std;
int main(void){
    int n;
    cin>>n;
    string ans;
    cin>>ans;
    int m;
    cin>>m;
    for(int i=0;i<m;i++){
        string str;
        cin>>str;
        int sum=0;
        for(int j=0;j<n;j++){
            if(str[j]==ans[j]){
                sum+=1;
            }
        }
        cout<<100/n*sum<<endl;
    }
}