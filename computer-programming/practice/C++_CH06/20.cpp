#include <bits\stdc++.h>
using namespace std;
int f(int m){
	if(m>0){
		return (2*f(m-1)-5);
	}else{
		return 3;
	}
}
int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int a;
			cin>>a;
			cout<<f(a)<<endl;
		}
	}
    return 0;
}

