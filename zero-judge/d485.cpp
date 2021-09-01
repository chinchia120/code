#include <iostream>
using namespace std;

int main() {
	int a,b,c;
	cin>>a>>b;
	c=b-a;
	if(a%2==0&&b%2==0){
		cout<<c/2+1<<endl;
	}else if(a%2!=0&&b%2!=0){
		cout<<c/2<<endl;
	}else{
		cout<<c/2+1<<endl;
	}
    return 0;
    }

