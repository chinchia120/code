#include <iostream>
using namespace std;

int main(void){
	char ptr[7][7]={};
	for(int i=1;i<6;i++){
		for(int j=1;j<6;j++){
			cin>>ptr[i][j];
			//cout<<ptr[i][j];
		}
		//cout<<endl;
	}
	for(int i=1;i<6;i++){
		for(int j=1;j<6;j++){
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
				cout<<cnt;
			}
		}
		cout<<endl;
	}
	return 0;
}
