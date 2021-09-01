#include <bits\stdc++.h>
using namespace std;
int mid1(int a,int b,int c){
	if(a>b && a>c){
		if(b>c){
			return b;
		}else{
			return c;
		}
	}else if(b>a && b>c){
		if(a>c){
			return a;
		}else{
			return c;
		}
	}else if(c>b && c>a){
		if(a>b){
			return a;
		}else{
			return b;
		}
	}
}
int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int a,b,c;
			cin>>a>>b>>c;
			cout<<mid1(a,b,c)<<endl;
		}
	}
    return 0;
}

