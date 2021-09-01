#include <iostream>
using namespace std;

int main() {
	int n;
	while(cin>>n){
		if(n==0){
			return 0;
		}else{
			for(int i=1;i<=n;i++){
				for(int j=n-i;j>=1;j--){
					cout<<"_";
				}
				for(int k=0;k<i;k++){
					cout<<"+";
				}
				cout<<endl;
			}
		}
	}
    return 0;
    }

