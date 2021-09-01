#include <bits\stdc++.h>
using namespace std;
void rf(int *a,int *b){
	int tmp=*b;
	for(int i=2;i<=tmp;i++){
		while(*a%i==0 && *b%i==0){
			*a=*a/i;
			*b=*b/i;
		}
	}
	cout<<*a<<"/"<<*b<<endl;
}
int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int a,b;
			cin>>a>>b;
			rf(&a,&b);
		}
	}
    return 0;
}

