#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            int m;
            cin>>m;
            int num[m];
            for(int j=0;j<m;j++){
                cin>>num[j];
            }
            for(int j=0;j<m;j++){
                for(int k=j+1;k<m;k++){
                    if(num[j]>num[k]){
                        int tmp=num[j];
                        num[j]=num[k];
                        num[k]=tmp;
                    }
                }
            }
            if(num[m-3]+num[m-2]>num[m-1] && num[m-3]+num[m-1]>num[m-2] && num[m-2]+num[m-1]>num[m-3]){
                cout<<"YES"<<endl;
            }else{
                cout<<"NO"<<endl;
            }
        }
    }
    return 0;
}