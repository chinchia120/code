#include <iostream>
using namespace std;

int main() {
	//question 1
	double a;
	int b;
	cout<<"input a positive float number A=";
	cin>>a;
	if(a<1){
		cout<<a<<endl;
	}else{
		b=a/1;
		cout<<a-b<<endl;
	}
    return 0;
}

