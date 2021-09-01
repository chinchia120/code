#include <iostream>
using namespace std;

int main() {
	int n,a;
	while(cin>>n){
		for(int i=0;i<n;i++){
			cin>>a;
			if(a==1){
				cout<<"No"<<endl;
				continue;
			} 
			int sum=0;
			for(int j=2;j<a;j++){
				if(a%j==0){
					sum=1;
				}
			}
			if(sum==1){
				cout<<"No"<<endl;
			}else if(sum==0){
				cout<<"Yes"<<endl;
			}
		}
	}
    return 0;
}

