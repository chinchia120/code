#include <iostream>
using namespace std;

int main() {
	int a,b,c,s,A;
	while(cin>>a>>b>>c){
		s=(a+b+c)/2;
		A=s*(s-a)*(s-b)*(s-c);
		cout<<A<<endl;
	}
    return 0;
    }

