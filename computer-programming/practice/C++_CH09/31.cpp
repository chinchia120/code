#include <bits\stdc++.h>
using namespace std;
void sumFactor(int m,int *sum){
	for(int i=1;i<m;i++){
		if(m%i==0){
			*sum=*sum+i;
		}
	}
}
int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int a,sum=0;
			cin>>a;
			sumFactor(a,&sum);
			cout<<sum<<endl;
		}
	}
    return 0;
}
