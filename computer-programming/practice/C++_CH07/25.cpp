#include <bits\stdc++.h>
using namespace std;
int f(int m){
	if(m>2){
		return (f(m-1)+f(m-2));
	}else{
		return 1;
	}
}
int main() {
	int n;
	cin>>n;
	for(int i=0;i<n;i++){
		int k;
		cin>>k;
		cout<<f(k)<<endl;
		
	}
    return 0;
}

