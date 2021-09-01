#include <iostream>
using namespace std;

int main() {
	int m,n;
	while(cin>>m>>n){
		if(m>n){
			cout<<n+1<<endl;
		}else if(m=n){
			cout<<m<<endl;
		}
	}
    return 0;
}

