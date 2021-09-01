#include <iostream>
using namespace std;

int main() {
	double n;
	int i=0;
	cin>>n;
	while(n>1){
		n=n/2;
		i=i+1;
	}
	cout<<i<<endl;
    return 0;
}

