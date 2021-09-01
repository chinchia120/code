#include <bits\stdc++.h>
using namespace std;
bool isPerfect(int m){
	int sum=0;
	for(int i=1;i<m;i++){
		if(m%i==0){
			sum=sum+i;
		}
	}
	if(sum==m){
		return true;
	}else{
		return false;
	}
}
int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int a,sum=0;
			cin>>a;
			if(isPerfect(a)){
				cout<<"Yes"<<endl;
			}else{
				cout<<"No"<<endl;
			}
			
		}
	}
    return 0;
}
