#include <bits\stdc++.h>
using namespace std;
bool isPrime(int m){
	int p=0;
	if(m==1){
		return false;
	}
	for(int i=2;i<m;i++){
		if(m%i==0){
			p=1;
		}
	}
	if(p==0){
		return true;
	}else{
		return false;
	}
}
int main() {
	int n;
	cin>>n;
	for(int i=0;i<n;i++){
		int k;
		cin>>k;
		if(isPrime(k)){
			cout<<"yes"<<endl;
		}else{
			cout<<"no"<<endl;
		}
	}
    return 0;
}

