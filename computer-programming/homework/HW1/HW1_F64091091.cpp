#include <iostream>
using namespace std;

int main() {
	int n,m;
	int b=0,s=0,sum=0,q=0;
	while(cin>>n>>m){
		if(m==0 && n==0){
			cout<<"end"<<endl;
			return 0;
		}else if(m<0 || n<0  || m>26 || n>26){
			cout<<"Please enter proper number."<<endl;
			continue;
		}else{
			//
			cout<<"A.  ";
			for(int i=0;i<n;i++){
				cout<<"#";
			}
			cout<<endl;
			//
			cout<<"B.  ";
			if(n<m){
				for(int i=n;i<m+1;i++){
					cout<<i<<" ";
				}
			}else if(n==m){
				cout<<n;
			}else{
				for(int i=m;i<n+1;i++){
					cout<<i<<" ";
				}
			}
			cout<<endl;
			//
			cout<<"C.  ";
			for(int i=0;i<m;i++){
				cout<<char(65+i);
			}
			cout<<endl;
			//
			cout<<"D.  ";
			for(int i=0;i<n;i++){
				cout<<char(122-i);
			}
			cout<<endl;
			//
			cout<<"E.  ";
			int l;
			if(m%2!=0){
				for(int i=1;i<m+1;i=i+2){
					cout<<i<<" ";
				}
			}else if(m%2==0){
				for(int i=m;i>1;i=i-2){
					cout<<i<<" ";
				}
			}	
			cout<<endl;
			//
			cout<<"F.  ";
			if(n<m){
				b=m;
				s=n;
			}else if(n==m){
				if(n%3==0){
					cout<<n;
				}
			}else if(m<n){
				b=n;
				s=m;
			}
			for(int i=s;i<b+1;i++){
				if(i%3==0){
					sum=sum+i;
				}
			}
			cout<<sum<<endl;
			//
			cout<<"G.  ";
			int g=m*n;
			for(int i=1;i<g+1;i++){
				if(g%i==0){
					cout<<i<<" ";
				}
			}
			cout<<endl;
			//
			cout<<"H.  ";
			int k=0;
			for(int i=1;i<101;i++){
				if(i%m==0 && i%n==0){
					cout<<i<<" ";
					k=1;
				}
			}
			if(k==0){
				cout<<"none";
			}
			cout<<endl;
			//
			cout<<"I.  ";
			if(n<m){
				b=m;
				s=n;
				while(b%s!=0){
					q=b%s;
					b=s;
					s=q;
				}
				while(b%s==0){
					cout<<s;
					break;
				}
			}else if(n==m){
				cout<<"1";
			}else if(m<n){
				b=n;
				s=m;
				while(b%s!=0){
					q=b%s;
					b=s;
					s=q;
				}
				while(b%s==0){
					cout<<s;
					break;
				}
			}
			cout<<endl;
			//
			cout<<"J.  ";
			for(int i=m;i<651;i++){
				if(i%m==0 && i%n==0){
					cout<<i;
					break;
				}
			}
			cout<<endl;	
			//
		}
	}
}

