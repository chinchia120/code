#include <iostream>
using namespace std;

int main(void){
	unsigned long long int a,b,c;
	while(cin>>a>>b>>c){
		if(a>b && a>c){
			if(b>c){
				if(a>b+c){
					cout<<"A"<<endl;
				}else if(a<b+c){
					cout<<"B"<<endl;
				}
			}else if(c>b){
				if(a>b+c){
					cout<<"A"<<endl;
				}else if(a<b+c){
					cout<<"C"<<endl;
				}
			}
		}else if(b>a && b>c){
			if(a>c){
				if(b>a+c){
					cout<<"B"<<endl;
				}else if(b<a+c){
					cout<<"A"<<endl;
				}
			}else if(c>a){
				if(b>a+c){
					cout<<"B"<<endl;
				}else if(b<a+c){
					cout<<"C"<<endl;
				}
			}
		}else if(c>a && c>b){
			if(a>b){
				if(c>a+b){
					cout<<"C"<<endl;
				}else if(c<a+b){
					cout<<"A"<<endl;
				}
			}else if(b>a){
				if(c>a+b){
					cout<<"C"<<endl;
				}else if(c<a+b){
					cout<<"B"<<endl;
				}
			}
		}
	}
	return 0;
}