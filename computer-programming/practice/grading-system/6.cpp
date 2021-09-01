#include <iostream>
using namespace std;

int main() {
	int n,a;
	while(cin>>n){
		for(int i=1;i<n+1;i++){
			cin>>a;
			cout<<"Case"<<i<<":";
			for(int j=1;j<a+1;j++){
				cout<<j;
			}
			cout<<endl;
		}
	}
    return 0;
}

