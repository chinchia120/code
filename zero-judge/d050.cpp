#include <iostream>
using namespace std;

int main() {
	int a;
	while(cin>>a){
		if(a>15){
			cout<<a-15<<endl;
		}else if(a==15){
			cout<<0<<endl;
		}else{
			cout<<a+24-15<<endl;
		}
	}
    return 0;
    }

