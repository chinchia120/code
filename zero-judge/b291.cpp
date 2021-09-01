#include <iostream>
using namespace std;

int main(){
    int n;
    string ans[4],guess[4];
    int A=0,B=0;        
    for(int i=0;i<4;i++){
        cin>>ans[i];
    }
    cin>>n;
    for(int i=0;i<n;i++){
        for(int j=0;j<4;j++){
            cin>>guess[j];
        }
        for(int l=0;l<4;l++){
            for(int k=0;k<4;k++){
                if(l==k && guess[l]==ans[k]){
                    A=A+1;
                }else if(l!=k && guess[l]==ans[k]){
                    B=B+1;
                }
            }
        }
        cout<<A<<"A"<<B<<"B"<<endl;
    }
    return 0;
}