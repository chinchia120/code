#include <bits\stdc++.h>
using namespace std;
char letter (int m){
	cout<<char(m+64)<<endl;
	return 1;
}
int main() {
	int n;
	cin>>n;
	for(int i=0;i<n;i++){
		int a;
		cin>>a;
		letter(a);
	}
    return 0;
}

