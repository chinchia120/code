#include <bits\stdc++.h>
using namespace std;
void gcd(int a,int b){
	int c;
	if(a>b){
		while(a%b!=0){
			c=a%b;
			a=b;
			b=c;
		}
		while(a%b==0){
			cout<<b<<endl;
			break;
		}
	}else{
		while(b%a!=0){
			c=b%a;
			b=a;
			a=c;
		}
		while(b%a==0){
			cout<<a<<endl;
			break;
		}
	}
} 
int main() {
	int n;
	cin>>n;
	for(int i=0;i<n;i++){
		int a,b;
		cin>>a>>b;
		gcd(a,b);
	}
    return 0;
}

