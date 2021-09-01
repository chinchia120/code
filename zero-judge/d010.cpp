#include <iostream>
using namespace std;

int main() {
	int n;
	while(cin>>n){
		int sum=0;
		for(int i=1;i<n;i++){
			if(n%i==0){
				sum=sum+i;
			}
		}
		if(sum>n){
			cout<<"盈數"<<endl;
		}else if(sum==n){
			cout<<"完全數"<<endl;
		}else if(sum<n){
			cout<<"虧數"<<endl;
		}
	}
    return 0;
}

