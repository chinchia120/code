#include <bits\stdc++.h>
using namespace std;
int digit(int m){
	int cnt=0;
	while(m!=0){
		m=m/10;
		cnt++;
	}
	while(m==0){
		return cnt;
	}
}
int main() {
	int n;
	cin>>n;
	for(int i=0;i<n;i++){
		int num;
		cin>>num;
		cout<<digit(num)<<endl;
	}
    return 0;
}

