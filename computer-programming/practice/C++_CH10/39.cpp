#include <iostream>
using namespace std;
int main(void){
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int a,b;
			char ch;
			cin>>a>>ch>>b;
			if(ch=='+'){
				cout<<a+b<<endl;
			}else if(ch=='-'){
				cout<<a-b<<endl;
			}else{
				cout<<a*b<<endl;
			}
		}
	}
	return 0;
}

