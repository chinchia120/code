#include <iostream>
using namespace std;

int main() {
	int m,n,a;
	cin>>m>>n;
	for(int i=0;i<m;i++){
		int sum=0;
		for(int j=0;j<n;j++){
			cin>>a;
			sum=sum+a;
		}
		cout<<sum<<endl;
	}
    return 0;
}
