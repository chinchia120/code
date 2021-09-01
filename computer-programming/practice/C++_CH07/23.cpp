#include <bits\stdc++.h>
using namespace std;
int lcm(int a,int b){
	int ans;
	for(int i=a;i<a*b;i++){
		if(i%a==0 && i%b==0){
			ans=i;
			return ans;
		}
	}
}
int main() {
	int n;
	cin>>n;
	for(int i=0;i<n;i++){
		int a,b;
		cin>>a>>b;
		cout<<lcm(a,b)<<endl;
	}
    return 0;
}

