#include <iostream>
using namespace std;

int main() {
	int n,a;
	while(cin>>n){
		for(int i=0;i<n;i++){
			cin>>a;
			int sum=0,b;
			while(a!=0){
				b=a%10;
				sum=sum+b;
				a=a/10;
			}
			while(a==0){
				cout<<sum<<endl;
				break;
			}
		}
	}
    return 0;
}

