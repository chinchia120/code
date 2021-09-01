#include <bits\stdc++.h>
using namespace std;
#define MAX(a,b,c) (a>b?a:b)>c?(a>b?a:b):c
int main() {
	int n;
	cin>>n;
	for(int i=0;i<n;i++){
		int a=0,b=0,c=0,max=0;
		cin>>a>>b>>c;
		max=MAX(a,b,c);
		cout<<max<<endl;
	}
    return 0;
}

