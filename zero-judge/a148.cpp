#include <iostream>
using namespace std;

int main() {
	int n;
	while(cin>>n){
		int s,t=0;
		for(int i=0;i<n;i++){
			cin>>s;
			t=t+s;
		}
		if(t<=59*n){
			cout<<"yes"<<endl;
		}else{
			cout<<"no"<<endl;
		}
	} 
    return 0;
    }

