#include <iostream>
using namespace std;

int main() {
	int n,sum=1;
	cin>>n;
	for(int i=1;i<n+1;i++){
		sum=sum*2;
	}
	cout<<sum<<endl;
    return 0;
}

