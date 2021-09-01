#include <iostream>
using namespace std;

int main() {
	long long int a,b;
	char ch;
	while(cin>>a>>ch>>b){
		if(ch=='+'){
			cout<<a+b<<endl;
		}else if(ch=='-'){
			cout<<a-b<<endl;
		}else if(ch=='*'){
			cout<<a*b<<endl;
		}else if(ch=='/'){
			cout<<a/b<<endl;
		}
	}
    return 0;
    }

