#include <iostream>
using namespace std;

int main(){
    string ans;
    while(cin>>ans){
        int n;
        cin>>n;
        for(int i=0;i<n;i++){
            string guess;
            cin>>guess;
            int A=0,B=0;
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    if(j==k && guess[j]==ans[k]){
                        A++;
                    }else if(j!=k && guess[j]==ans[k]){
                        B++;
                    }
                }
            }
            cout<<A<<"A"<<B<<"B"<<endl; 
        }
    }
    return 0;
}