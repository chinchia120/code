#include <iostream>
using namespace std;

int main(void){
    int n;
    cin>>n;
	char ptr[n+2][n+2]={};
	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			cin>>ptr[i][j];
		}
	}
	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			if(ptr[i][j]=='*'){
				cout<<"*";
			}else{
				int cnt=0;
				for(int a=i-1;a<i+2;a++){
					for(int b=j-1;b<j+2;b++){
						if(ptr[a][b]=='*'){
							cnt++;
						}
					}
				}
                if(cnt==0){
                    cout<<"-";
                }else{
                    cout<<cnt;
                }
			}
		}
		cout<<endl;
	}
	return 0;
}
