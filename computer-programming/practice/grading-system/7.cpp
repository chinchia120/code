#include <iostream>
using namespace std;

int main() {
	int n,a;
	while(cin>>n){
		for(int i=0;i<n;i++){
			cin>>a;
			if(a%2==0){
				cout<<-a/2<<endl;
			}else if(a%2!=0){
				cout<<(a+1)/2<<endl;
			}
		}
	}
    return 0;
}

